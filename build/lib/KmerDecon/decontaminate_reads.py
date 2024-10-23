# src/KmerDecon/decontaminate_reads.py

import argparse
import os
import sys
from KmerDecon.bloom_filter import BloomFilter
from KmerDecon.utils import generate_kmers
from Bio import SeqIO
import csv

def main():
    parser = argparse.ArgumentParser(
        description='Decontaminate sequencing reads using Bloom filters.'
    )
    parser.add_argument('-i', '--input-reads', required=True,
                        help='Input FASTQ file or directory containing FASTQ files.')
    parser.add_argument('-b', '--bloom-filter', required=True,
                        help='Bloom filter file or directory containing bloom filters for contamination.')
    parser.add_argument('-t', '--threshold', type=float, default=0.5,
                        help='Contamination threshold (default: 0.5).')
    parser.add_argument('-k', '--kmer-length', type=int,
                        help='Length of k-mers. If not provided, it will use the k-mer length from the Bloom filters.')
    parser.add_argument('-o', '--output-dir', required=True,
                        help='Output directory for decontaminated reads or states.csv.')
    parser.add_argument('--mode', choices=['filter', 'states'], default='filter',
                        help='Mode of operation: "filter" to filter reads, "states" to generate a states.csv file.')

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    elif not os.path.isdir(args.output_dir):
        print(f"Error: {args.output_dir} is not a directory.")
        sys.exit(1)

    print("Loading Bloom filters...")
    bloom_filters = []

    if os.path.isfile(args.bloom_filter):
        bloom_filter = BloomFilter.load(args.bloom_filter)
        filter_name = os.path.basename(args.bloom_filter).split('.')[0]
        bloom_filters.append((filter_name, bloom_filter))
    elif os.path.isdir(args.bloom_filter):
        for filename in os.listdir(args.bloom_filter):
            if filename.endswith('.bf'):  # Assuming bloom filter files have .bf extension
                bf_path = os.path.join(args.bloom_filter, filename)
                bloom_filter = BloomFilter.load(bf_path)
                filter_name = os.path.basename(filename).split('.')[0]
                bloom_filters.append((filter_name, bloom_filter))
    else:
        print(f"Error: {args.bloom_filter} is not a valid file or directory.")
        sys.exit(1)

    print(f"Loaded {len(bloom_filters)} bloom filter(s).")

    kmer_lengths = set()

    for filter_name, bloom_filter in bloom_filters:
        kmer_lengths.add(bloom_filter.kmer_length)

    if args.kmer_length:
        k = args.kmer_length
        if len(kmer_lengths) > 1 or (len(kmer_lengths) == 1 and k not in kmer_lengths):
            print(f"Warning: Provided k-mer length ({k}) does not match the k-mer length(s) used in the Bloom filters ({kmer_lengths}).")
    else:
        if len(kmer_lengths) == 1:
            k = kmer_lengths.pop()
            print(f"Using k-mer length {k} from the Bloom filters.")
        else:
            print(f"Error: Multiple k-mer lengths found in Bloom filters ({kmer_lengths}). Please specify --kmer-length.")
            sys.exit(1)

    input_files = []
    if os.path.isfile(args.input_reads):
        input_files.append(args.input_reads)
    elif os.path.isdir(args.input_reads):
        for filename in os.listdir(args.input_reads):
            if filename.endswith('.fastq') or filename.endswith('.fq'):
                input_files.append(os.path.join(args.input_reads, filename))
    else:
        print(f"Error: {args.input_reads} is not a valid file or directory.")
        sys.exit(1)

    if args.mode == 'states':
        states_filename = os.path.join(args.output_dir, 'states.csv')
        csv_file = open(states_filename, 'w', newline='')
        fieldnames = ['input_file']
        for filter_name, _ in bloom_filters:
            fieldnames.append(f"{filter_name}_avgFrac")
            fieldnames.append(f"{filter_name}_percentPassing")
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

    for input_file in input_files:
        print(f"Processing input file: {input_file}")
        if args.mode == 'filter':
            output_filename = os.path.join(args.output_dir, os.path.basename(input_file))
            with open(output_filename, 'w') as out_f:
                total_reads = 0
                kept_reads = 0
                for record in SeqIO.parse(input_file, "fastq"):
                    total_reads += 1
                    seq = str(record.seq).upper()
                    if len(seq) < k:
                        continue
                    num_kmers = len(seq) - k + 1
                    if num_kmers == 0:
                        continue
                    kmers = list(generate_kmers(seq, k))
                    keep_read = True
                    for filter_name, bloom_filter in bloom_filters:
                        count = sum(1 for kmer in kmers if kmer in bloom_filter)
                        fraction = count / num_kmers
                        if fraction >= args.threshold:
                            keep_read = False
                            break
                    if keep_read:
                        SeqIO.write(record, out_f, "fastq")
                        kept_reads += 1
                print(f"File {input_file}: {kept_reads}/{total_reads} reads kept. Output written to {output_filename}")
        elif args.mode == 'states':
            total_reads = 0
            fractions_sum = {filter_name: 0.0 for filter_name, _ in bloom_filters}
            passing_counts = {filter_name: 0 for filter_name, _ in bloom_filters}
            for record in SeqIO.parse(input_file, "fastq"):
                total_reads += 1
                seq = str(record.seq).upper()
                if len(seq) < k:
                    continue
                num_kmers = len(seq) - k + 1
                if num_kmers == 0:
                    continue 
                kmers = list(generate_kmers(seq, k))
                for filter_name, bloom_filter in bloom_filters:
                    count = sum(1 for kmer in kmers if kmer in bloom_filter)
                    fraction = count / num_kmers
                    fractions_sum[filter_name] += fraction
                    if fraction < args.threshold:
                        passing_counts[filter_name] += 1
            row = {'input_file': os.path.basename(input_file)}
            for filter_name, _ in bloom_filters:
                avg_fraction = fractions_sum[filter_name] / total_reads if total_reads > 0 else 0
                percent_passing = (passing_counts[filter_name] / total_reads * 100) if total_reads > 0 else 0
                row[f"{filter_name}_avgSimilarity"] = avg_fraction
                row[f"{filter_name}_percentReadsPassing"] = percent_passing
            writer.writerow(row)
            print(f"File {input_file}: Processed {total_reads} reads.")
    if args.mode == 'states':
        csv_file.close()
        print(f"States written to {states_filename}")

if __name__ == "__main__":
    main()