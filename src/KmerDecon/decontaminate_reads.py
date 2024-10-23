# src/KmerDecon/decontaminate_reads.py

import argparse
from KmerDecon.bloom_filter import BloomFilter
from KmerDecon.utils import generate_kmers
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(
        description='Decontaminate sequencing reads using a Bloom filter.'
    )
    parser.add_argument('-i', '--input-reads', required=True,
                        help='Input FASTQ file with sequencing reads.')
    parser.add_argument('-b', '--bloom-filter', required=True,
                        help='Bloom filter file for contamination.')
    parser.add_argument('-t', '--threshold', type=float, default=0.5,
                        help='Contamination threshold (default: 0.5).')
    parser.add_argument('-k', '--kmer-length', type=int, required=True,
                        help='Length of k-mers.')
    parser.add_argument('-o', '--output-reads', required=True,
                        help='Output FASTQ file for decontaminated reads.')

    args = parser.parse_args()

    # Load the Bloom filter
    print("Loading Bloom filter...")
    bloom_filter = BloomFilter.load(args.bloom_filter)
    print("Bloom filter loaded.")

    with open(args.output_reads, 'w') as out_f:
        print("Starting decontamination...")
        total_reads = 0
        kept_reads = 0
        for record in SeqIO.parse(args.input_reads, "fastq"):
            total_reads += 1
            seq = str(record.seq).upper()
            kmers = list(generate_kmers(seq, args.kmer_length))
            if not kmers:
                continue  # Skip reads shorter than k
            count = sum(1 for kmer in kmers if kmer in bloom_filter)
            fraction = count / len(kmers)
            if fraction < args.threshold:
                # Keep the read
                SeqIO.write(record, out_f, "fastq")
                kept_reads += 1
            # Else, the read is discarded (contaminated)
        print(f"Decontamination completed. {kept_reads}/{total_reads} reads kept.")

if __name__ == "__main__":
    main()