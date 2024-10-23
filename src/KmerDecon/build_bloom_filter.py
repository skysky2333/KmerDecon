# src/KmerDecon/build_bloom_filter.py

import argparse
from KmerDecon.bloom_filter import BloomFilter
from KmerDecon.utils import generate_kmers
from Bio import SeqIO
from hyperloglog import HyperLogLog
import math
import statistics

def estimate_unique_kmers(contamination_fasta: str, k: int) -> int:
    """
    Estimate the number of unique k-mers in the contamination sequences using HyperLogLog.

    Args:
        contamination_fasta (str): Path to the contamination FASTA file.
        k (int): Length of k-mers.

    Returns:
        int: Estimated number of unique k-mers.
    """
    print(f"Estimating the number of unique {k}-mers in contamination sequences using HyperLogLog...")
    hll = HyperLogLog(0.01)  # 1% relative error
    for record in SeqIO.parse(contamination_fasta, "fasta"):
        seq = str(record.seq).upper()
        for kmer in generate_kmers(seq, k):
            hll.add(kmer)
    n_unique = int(len(hll))
    print(f"Estimated {n_unique} unique {k}-mers.")
    return n_unique

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
    for record in SeqIO.parse(contamination_fasta, "fasta"):
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
    # Ensure k is odd (common practice)
    if suggested_k % 2 == 0:
        suggested_k += 1
    print(f"Suggested k-mer length: {suggested_k}")
    return suggested_k

def main():
    parser = argparse.ArgumentParser(
        description='Build a Bloom filter from contamination sources.'
    )
    parser.add_argument('-c', '--contamination-fasta', required=True,
                        help='FASTA file with contamination sequences.')
    parser.add_argument('-k', '--kmer-length', type=int,
                        help='Length of k-mers. If not provided, it will be determined automatically.')
    parser.add_argument('-o', '--output-filter', required=True,
                        help='Output Bloom filter file.')
    parser.add_argument('-p', '--false-positive-rate', type=float, default=0.001,
                        help='Desired false positive rate (default: 0.001).')
    parser.add_argument('-e', '--expected-elements', type=int,
                        help='Expected number of unique k-mers. If not provided, it will be estimated.')
    parser.add_argument('-m', '--max-memory', type=float,
                        help='Maximum memory in GB for the Bloom filter. Overrides false positive rate if set.')

    args = parser.parse_args()

    # Determine k-mer length if not provided
    if args.kmer_length:
        k = args.kmer_length
    else:
        k = determine_best_kmer_length(args.contamination_fasta)

    if args.expected_elements:
        n_unique = args.expected_elements
    else:
        n_unique = estimate_unique_kmers(args.contamination_fasta, k)

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
    print(f"Bloom filter size: {bloom_size_bytes / (1024 ** 3):.4f} GB")
    print(f"Number of hash functions: {bloom_filter.hash_count}")

    print("Building Bloom filter...")
    for record in SeqIO.parse(args.contamination_fasta, "fasta"):
        seq = str(record.seq).upper()
        for kmer in generate_kmers(seq, k):
            bloom_filter.add(kmer)
    print("Bloom filter built.")

    bloom_filter.save(args.output_filter)
    print(f"Bloom filter saved to {args.output_filter} and {args.output_filter}.params")

if __name__ == "__main__":
    main()