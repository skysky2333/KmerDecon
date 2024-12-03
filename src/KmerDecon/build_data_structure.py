# src/KmerDecon/build_bloom_filter.py
import argparse
from KmerDecon.bloom_filter import BloomFilter
from cms import CountMinSketch
from KmerDecon.utils import generate_kmers
from Bio import SeqIO
from hyperloglog import HyperLogLog
import math
import statistics
from tqdm import tqdm

def estimate_unique_kmers(contamination_fasta: str, k: int, exclude_filter: BloomFilter = None) -> int:
    """
    Estimate the number of unique k-mers in the contamination sequences using HyperLogLog.

    Args:
        contamination_fasta (str): Path to the contamination FASTA file.
        k (int): Length of k-mers.
        exclude_filter (BloomFilter): excluded filter used.

    Returns:
        int: Estimated number of unique k-mers.
    """
    print(f"Estimating the number of unique {k}-mers in contamination sequences using HyperLogLog...")
    hll = HyperLogLog(0.01)  # 1% relative error
    total_kmers = 0
    for record in tqdm(SeqIO.parse(contamination_fasta, "fasta"), desc="Estimating the number of unique k-mers"):
        seq = str(record.seq).upper()
        for kmer in generate_kmers(seq, k):
            total_kmers += 1
            if exclude_filter and kmer in exclude_filter:
                continue
            hll.add(kmer)
    n_unique = int(len(hll))
    print(f"Estimated {n_unique} unique {k}-mers out of {total_kmers} total k-mers.")
    return n_unique, total_kmers

def determine_best_kmer_length(contamination_fasta: str) -> int:
    """
    Automatically determine the best k-mer length based on sequence lengths.

    Args:
        contamination_fasta (str): Path to the contamination FASTA file.

    Returns:
        int: Suggested k-mer length.
    """
    print("Determining the best k-mer length...")
    lengths = []
    for record in tqdm(SeqIO.parse(contamination_fasta, "fasta"),desc="Determining the best k-mer length"):
        lengths.append(len(record.seq))
    if not lengths:
        print("No sequences found in the contamination FASTA file.")
        exit(1)
    median_length = int(statistics.median(lengths))
    suggested_k = median_length // 2
    # Enforce reasonable bounds
    if suggested_k < 21:
        suggested_k = 21
    elif suggested_k > 127:
        suggested_k = 127
    # Ensure k is odd
    if suggested_k % 2 == 0:
        suggested_k += 1
    print(f"Suggested k-mer length: {suggested_k}")
    return suggested_k

def main():
    parser = argparse.ArgumentParser(description="Build or load data structures (Bloom Filter or CMS) from contamination sequences.")
    parser.add_argument('-c', '--contamination-fasta', required=True, 
                        help='FASTA file with contamination sequences.')
    parser.add_argument('-k', '--kmer-length', type=int, 
                        help='Length of k-mers. If not provided, it will be determined automatically.')
    parser.add_argument('-o', '--output-filter', required=True, 
                        help='Output file for the data structure (either Bloom filter or CMS).')
    parser.add_argument('-s', '--data-structure', choices=['bloom', 'cms'], required=True,
                        help='Choose whether to build a Bloom filter or CountMinSketch.')
    parser.add_argument('-p', '--false-positive-rate', type=float, default=0.01, 
                        help='Desired false positive rate for Bloom filter (default: 0.01).')
    parser.add_argument('-e', '--expected-elements', type=int, 
                        help='Expected number of unique k-mers. If not provided, it will be estimated.')
    parser.add_argument('-m', '--max-memory', type=float,
                        help='Maximum memory in GB for the Bloom filter. Overrides false positive rate if set.')
    parser.add_argument('-x', '--exclude-filter', 
                        help='Bloom filter or CountMinSketch file to exclude kmers from.')
    parser.add_argument('-w', '--width-of-CountMinSketch', type=int,
                        help='The width of Count_mint_sketch')
    parser.add_argument('-d', '--depth-of-CountMinSketch', type=int,
                        help='The depth of Count_mint_sketch')
    parser.add_argument('-r', '--error-rate', type=float,
                        help='The error rate of Count-Min Sketch (default: 0.01)')


    args = parser.parse_args()
    # Determine k-mer length if not provided
    if args.data_structure=='bloom':
        if args.exclude_filter:
            print("Loading exclude bloom filter...")
            exclude_filter = BloomFilter.load(args.exclude_filter)
            k = exclude_filter.kmer_length
            print(f"Using k-mer length {k} from the exclude bloom filter.")
        else:
            if args.kmer_length:
                k = args.kmer_length
            else:
                k = determine_best_kmer_length(args.contamination_fasta)

        if args.expected_elements:
            n_unique = args.expected_elements
            total_kmers = None
        else:
            n_unique, total_kmers = estimate_unique_kmers(args.contamination_fasta, k, exclude_filter if args.exclude_filter else None)
            
        if args.max_memory:
            # Calculate false positive rate based on max memory
            max_bits = args.max_memory * 8 * (1024 ** 3)  # Convert GB to bits
            p = math.exp(- (max_bits * (math.log(2) ** 2)) / n_unique)
            false_positive_rate = p
            print(f"Adjusted false positive rate to {false_positive_rate:.6f} based on max memory {args.max_memory} GB.")
        else:
            false_positive_rate = args.false_positive_rate

        bloom_filter = BloomFilter(n_unique, false_positive_rate, k)

        bloom_size_bytes = bloom_filter.size / 8
        print(f"Bloom filter size: {bloom_size_bytes / (1024 ** 3):.4f} GB, est. file size {bloom_size_bytes / (1024 ** 3)*30:.4f} MB")
        print(f"Number of hash functions: {bloom_filter.hash_count}")

        print("Building Bloom filter...")
        total_kmers = 0
        unique_kmers = 0
        for record in tqdm(SeqIO.parse(args.contamination_fasta, "fasta"), desc="Building Bloom filter"):
            seq = str(record.seq).upper()
            for kmer in generate_kmers(seq, k):
                total_kmers += 1
                if args.exclude_filter and kmer in exclude_filter:
                    continue
                bloom_filter.add(kmer)
                unique_kmers += 1
        if total_kmers > 0:
            percent_unique = (unique_kmers / total_kmers) * 100
            print(f"{percent_unique:.2f}% of k-mers are unique and encoded in the Bloom filter.")
        else:
            print("No k-mers were processed.")

        bloom_filter.save(args.output_filter)
        print(f"Bloom filter saved to {args.output_filter} and {args.output_filter}.params")



    elif args.data_structure=='cms':
        if args.exclude_filter:
            print("Loading exclude CountMinSketch...")
            exclude_filter = CountMinSketch.load(args.exclude_filter)
            k = exclude_filter.kmer_length
            print(f"Using k-mer length {k} from the exclude CountMinSketch.")
        else:
            if args.kmer_length:
                k = args.kmer_length
            else:
                k = determine_best_kmer_length(args.contamination_fasta)
        if args.expected_elements:
            n_unique = args.expected_elements
            total_kmers = None
        else:
            n_unique, total_kmers = estimate_unique_kmers(args.contamination_fasta, k, exclude_filter if args.exclude_filter else None)
        

        if args.width_of_CountMinSketch and args.depth_of_CountMinSketch:
            w=args.width_of_CountMinSketch
            d=args.depth_of_CountMinSketch
        elif args.error_rate:
            error_rate=args.error_rate
            w,d=CountMinSketch._optimal_params(n_unique,error_rate)
        elif args.error_rate== None:
            w,d=CountMinSketch._optimal_params(n_unique)
        cms = CountMinSketch(w,d,k)
        cms_size_bytes = cms.size / 8
        print(f" size: {cms_size_bytes / (1024 ** 3):.4f} GB, est. file size {cms_size_bytes / (1024 ** 3)*30:.4f} MB")
        print("Building Count_Min_Sketch...")

        total_kmers = 0
        unique_kmers = 0
        for record in tqdm(SeqIO.parse(args.contamination_fasta, "fasta"), desc="Building Count_Mint_Sketch"):
            seq = str(record.seq).upper()
            for kmer in generate_kmers(seq, k):
                total_kmers += 1
                if args.exclude_filter and kmer in exclude_filter:
                    continue
                cms.add(kmer)
                unique_kmers += 1
        if total_kmers > 0:
            percent_unique = (unique_kmers / total_kmers) * 100
            print(f"{percent_unique:.2f}% of k-mers are unique and encoded in the Count_Mint_Sketch.")
        else:
            print("No k-mers were processed.")

        cms.save(args.output_filter)
        print(f"Count_Mint_Sketch saved to {args.output_filter} and {args.output_filter}.params")
        a=cms.count("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
        print(a)

if __name__ == "__main__":
    main()