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
    parser.add_argument('-k', '--kmer-length', type=int,
                        help='Length of k-mers. If not provided, it will use the k-mer length from the Bloom filter.')
    parser.add_argument('-o', '--output-reads', required=True,
                        help='Output FASTQ file for decontaminated reads.')

    args = parser.parse_args()

    print("Loading Bloom filter...")
    bloom_filter = BloomFilter.load(args.bloom_filter)
    print("Bloom filter loaded.")

    # Determine k-mer length
    if args.kmer_length:
        k = args.kmer_length
        if k != bloom_filter.kmer_length:
            print(f"Warning: Provided k-mer length ({k}) does not match the k-mer length used in the Bloom filter ({bloom_filter.kmer_length}).")
    else:
        k = bloom_filter.kmer_length
        print(f"Using k-mer length {k} from the Bloom filter.")

    with open(args.output_reads, 'w') as out_f:
        print("Starting decontamination...")
        total_reads = 0
        kept_reads = 0
        for record in SeqIO.parse(args.input_reads, "fastq"):
            total_reads += 1
            seq = str(record.seq).upper()
            if len(seq) < k:
                continue  # Skip reads shorter than k
            kmers = generate_kmers(seq, k)
            count = sum(1 for kmer in kmers if kmer in bloom_filter)
            num_kmers = len(seq) - k + 1
            fraction = count / num_kmers
            if fraction < args.threshold:
                SeqIO.write(record, out_f, "fastq")
                kept_reads += 1
        print(f"Decontamination completed. {kept_reads}/{total_reads} reads kept.")

if __name__ == "__main__":
    main()