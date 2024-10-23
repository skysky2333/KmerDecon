# src/KmerDecon/build_bloom_filter.py

import argparse
from KmerDecon.bloom_filter import BloomFilter
from KmerDecon.utils import generate_kmers
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(
        description='Build a Bloom filter from contamination sources.'
    )
    parser.add_argument('-c', '--contamination-fasta', required=True,
                        help='FASTA file with contamination sequences.')
    parser.add_argument('-k', '--kmer-length', type=int, required=True,
                        help='Length of k-mers.')
    parser.add_argument('-o', '--output-filter', required=True,
                        help='Output Bloom filter file.')
    parser.add_argument('-e', '--expected-elements', type=int, default=int(1e8),
                        help='Expected number of k-mers (default: 1e8).')
    parser.add_argument('-f', '--false-positive-rate', type=float, default=0.001,
                        help='Desired false positive rate (default: 0.001).')

    args = parser.parse_args()

    bloom_filter = BloomFilter(args.expected_elements, args.false_positive_rate)

    print("Building Bloom filter...")
    for record in SeqIO.parse(args.contamination_fasta, "fasta"):
        seq = str(record.seq).upper()
        for kmer in generate_kmers(seq, args.kmer_length):
            bloom_filter.add(kmer)
    print("Bloom filter built.")

    bloom_filter.save(args.output_filter)
    print(f"Bloom filter saved to {args.output_filter} and {args.output_filter}.params")

if __name__ == "__main__":
    main()