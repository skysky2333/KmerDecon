# src/KmerDecon/build_bloom_filter.py
import argparse
from KmerDecon.bloom_filter import BloomFilter
from KmerDecon.cuckoofilter import CuckooFilter
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

def main():
    parser = argparse.ArgumentParser(description="Build or load data structures (Bloom Filter or CMS) from contamination sequences.")
    parser.add_argument('-c', '--contamination-fasta', required=True, 
                        help='FASTA file with contamination sequences.')
    parser.add_argument('-k', '--kmer-length', type=int, default=31,
                        help='Length of k-mers. Default is 31.')
    parser.add_argument('-o', '--output-filter', required=True, 
                        help='Output file for the data structure (either Bloom filter or CMS).')
    parser.add_argument('-s', '--data-structure', choices=['bloom', 'cuckoo'], required=True,
                        help='Choose whether to build a Bloom filter or Cuckoo filter.')
    parser.add_argument('-p', '--false-positive-rate', type=float, default=0.01, 
                        help='Desired false positive rate for Bloom filter (default: 0.01).')
    parser.add_argument('-e', '--expected-elements', type=int, 
                        help='Expected number of unique k-mers. If not provided, it will be estimated.')
    parser.add_argument('-m', '--max-memory', type=float,
                        help='Maximum memory in GB for the Bloom filter. Overrides false positive rate if set.')
    parser.add_argument('-x', '--exclude-filter', 
                        help='Bloom filter or Cuckoo filter file to exclude kmers from.')
    parser.add_argument('-cap', '--capacity-of-cuckoofilter', type=int,
                        help='The capacity of cuckoo filter')
    parser.add_argument('-f', '--fingerprint-size-of-cuckoofilter', type=int,
                        help='The fingerprint size of cuckoo filter')

    args = parser.parse_args()
    
    '''if choose to build bloom filter, 
       load exclude bloom filter, 
       choose appropriate kmer length (dafault 31) and false positive, 
       build new bloom filter
    '''
    
    if args.data_structure=='bloom':
        if args.exclude_filter:
            print("Loading exclude bloom filter...")
            exclude_filter = BloomFilter.load(args.exclude_filter)
            k = exclude_filter.kmer_length
            print(f"Using k-mer length {k} from the exclude bloom filter.")
        else:
            k = args.kmer_length

        if args.expected_elements:
            n_unique = args.expected_elements
            total_kmers = None
        else:
            #estimate number of unique kmer and total kmer
            n_unique, total_kmers = estimate_unique_kmers(args.contamination_fasta, k, exclude_filter if args.exclude_filter else None)
            
        if args.max_memory:
            # Calculate false positive rate based on max memory
            max_bits = args.max_memory * 8 * (1024 ** 3)  # Convert GB to bits
            p = math.exp(- (max_bits * (math.log(2) ** 2)) / n_unique)
            false_positive_rate = p
            print(f"Adjusted false positive rate to {false_positive_rate:.6f} based on max memory {args.max_memory} GB.")
        else:
            false_positive_rate = args.false_positive_rate

        #build new bloom filter
        bloom_filter = BloomFilter(n_unique, false_positive_rate, k)

        # calculate the size of the Bloom filter
        bloom_size_bytes = bloom_filter.size / 8
        print(f"Bloom filter size: {bloom_size_bytes / (1024 ** 3):.4f} GB")
        print(f"Number of hash functions: {bloom_filter.hash_count}")

        print("Building Bloom filter...")
        total_kmers = 0
        unique_kmers = 0
        #use SeqIO get each sequence from fasta file
        for record in tqdm(SeqIO.parse(args.contamination_fasta, "fasta"), desc="Building Bloom filter"):
            seq = str(record.seq).upper()
            for kmer in generate_kmers(seq, k):
                total_kmers += 1
                if args.exclude_filter and kmer in exclude_filter:
                    continue
                bloom_filter.add(kmer)#add kmer to bloom filter
                unique_kmers += 1 #count the number of unique kmers
        if total_kmers > 0:
            percent_unique = (unique_kmers / total_kmers) * 100
            print(f"{percent_unique:.2f}% of k-mers are unique and encoded in the Bloom filter.")
        else:
            print("No k-mers were processed.")

        bloom_filter.save(args.output_filter)#save bloom filter as output file
        print(f"Bloom filter saved to {args.output_filter} and {args.output_filter}.params")


        '''
        if choose build cuckoo filter
        load exclude cukoo filter
        get appropriate kmer length (default 31) and false positive
        build new cuckoo filter

        '''
    elif args.data_structure=='cuckoo':
        if args.exclude_filter:
            print("Loading exclude cuckoo filter...")
            exclude_filter = CuckooFilter.load(args.exclude_filter)
            k = exclude_filter.kmer_length
            print(f"Using k-mer length {k} from the exclude CuckooFilter.")
        else:
            if args.kmer_length:
                k = args.kmer_length
            else:
                k = 31 #set default kmer length as 31

        if args.expected_elements:
            n_unique = args.expected_elements
            total_kmers = None
        else:
            #estimate the number of unique kmer and total kmer
            n_unique, total_kmers = estimate_unique_kmers(args.contamination_fasta, k, exclude_filter if args.exclude_filter else None)

        '''calculate the bucket size and estimate the size of fingerprint and capacity
            the euqation get from https://stackoverflow.com/questions/57555236/how-to-size-a-cuckoo-filter
        '''
        if args.false_positive_rate<0.002:
            bucket_size=4
        else:
            bucket_size=2
        
        print(f"Using bucket size of {bucket_size}")

        if args.capacity_of_cuckoofilter and args.fingerprint_size_of_cuckoofilter:
            capacity=args.capacity_of_cuckoofilter
            fingerprint_size=args.fingerprint_size_of_cuckoofilter
        else:
            fingerprint_size=int(math.log2(1 / args.false_positive_rate)+math .log2(2*bucket_size))
            if bucket_size==4:
                capacity=int(total_kmers/0.95)
            else:
                capacity=int(total_kmers/0.84)
        #build new cuckoo filter
        cuckoo = CuckooFilter(capacity,fingerprint_size,k, bucket_size)
        cuckoo_size_bytes = cuckoo.__sizeof__()
        print(f" size: {cuckoo_size_bytes / (1024 ** 3):.4f} GB, est. file size {cuckoo_size_bytes / (1024 ** 2):.4f} MB")
        print("Building Cuckoo Filter...")

        total_kmers = 0
        unique_kmers = 0
        #use SeqIO get the sequence from fasta file
        for record in tqdm(SeqIO.parse(args.contamination_fasta, "fasta"), desc="Building Cuckoo filter"):
            seq = str(record.seq).upper()
            for kmer in generate_kmers(seq, k):
                total_kmers += 1#count the number of total kmer
                if args.exclude_filter and kmer in exclude_filter:
                    continue
                cuckoo.insert(kmer)#add kmer in cuckoo filter
                unique_kmers += 1#count the number of unique kmer
        if total_kmers > 0:
            percent_unique = (unique_kmers / total_kmers) * 100
            print(f"{percent_unique:.2f}% of k-mers are unique and encoded in the cuckoo filter.")
        else:
            print("No k-mers were processed.")
        cuckoo.save(args.output_filter)#save cuckoo filter as output file
        print(f"cuckoo filter saved to {args.output_filter} and {args.output_filter}.params")

if __name__ == "__main__":
    main()
