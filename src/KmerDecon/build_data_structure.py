# src/KmerDecon/build_bloom_filter.py
import argparse
from KmerDecon.bloom_filter import BloomFilter
from KmerDecon.cuckoofilter import BCuckooFilter
from KmerDecon.utils import generate_kmers
from Bio import SeqIO
from hyperloglog import HyperLogLog
import math
from tqdm import tqdm
from multiprocessing import Pool
import copy
from functools import partial

# Global variables used by workers in parrallel process
exclude_filter = None
bloom_filter = None

def init_hll_worker(global_exclude_filter):
    """
    Initializes the HyperLogLog worker for parallel proccessing.
    Each worker gets its own copy of the exclude filter.
    """
    global exclude_filter
    exclude_filter = copy.deepcopy(global_exclude_filter)

def init_bf_worker(global_exclude_filter, n_unique, false_positive_rate, k):
    """
    Initializes the Bloom filter worker for parallel proccessing.
    Each worker gets its own copy of the exclude filter.
    Each worker creates its own Bloom filter.
    """
    global exclude_filter
    exclude_filter = copy.deepcopy(global_exclude_filter)
    global bloom_filter
    bloom_filter = BloomFilter(n_unique, false_positive_rate, k)


def process_seq_for_kmers(seq, k):
    """
    This is only used when using multiple cores.
    Estimate the number of unique k-mers in the sequence using HyperLogLog.

    Args:
        seq SeqRecord: Sequence to calculate for.
        k (int): Length of k-mers.

    Returns:
        int: The number of kmers in the sequence.
        HyperLogLog: The HyperLogLog generated from this sequence.
    """
    local_hll = HyperLogLog(0.01)
    total_kmers = 0
    seq = str(seq.seq).upper()
    for kmer in generate_kmers(seq, k):
        total_kmers += 1
        if exclude_filter and kmer in exclude_filter:
            continue
        local_hll.add(kmer)  # Add k-mer to the local HyperLogLog instance
    return total_kmers, local_hll

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
    return n_unique

def estimate_unique_kmers_parallel(records, k, n, exclude_filter_in=None):
    """
    Estimate the number of unique k-mers in the contamination sequences using HyperLogLog with multiprocessing.

    Args:
        records List[SeqRecord]: The list of sequences.
        k (int): Length of k-mers.
        exclude_filter (BloomFilter): excluded filter used.

    Returns:
        int: Estimated number of unique k-mers.
    """
    
    print(f"Got {len(records)} sequences.")

    total_kmers = 0

    global_hll = HyperLogLog(0.01)  

    total_seq = 0

    with Pool(processes=n, initializer=init_hll_worker, initargs=(exclude_filter_in,)) as pool:
        with tqdm(total=len(records), desc="Estimating unique k-mers") as pbar:
             func = partial(process_seq_for_kmers, k=k)
             for result in pool.imap_unordered(func, records):
                total_kmers += result[0]
                # merge the local hll into global hlll
                global_hll.update(result[1])
                total_seq += 1
                pbar.update(1)

    n_unique = int(len(global_hll)) 

    print(f"Estimated unique k-mers: {n_unique}")
    print(f"Total k-mers: {total_kmers}")
    print(f"Total sequences: {total_seq}")
    
    return n_unique

def process_sequence_chunk(chunk, k):
    """
    Used in parallel processing.
    Process each sequence record to add k-mers to the Bloom filter.

    Args:
        chunk (List[SeqRecord]): The list of sequences to build the Bloom filter for.
        k (int): Length of k-mers.

    Returns:
        total_kmers (int): Total kmers in the input chunk.
        bloom_filter (BloomFilter): The Bloom filter constructed by the worker.
        total_seq (int): The number of sequences.
        kmers_added (int): The number of kmers not excluded by the exclude filter.
    """
    total_kmers = 0
    total_seq = 0
    kmers_added = 0

    for seq in tqdm(chunk, desc="Adding k-mers from chunk"):
        seq = str(seq.seq).upper()
        
        for kmer in generate_kmers(seq, k):
            total_kmers += 1
            # if kmer in unique_kmers_set:
            #     continue
            if exclude_filter and kmer in exclude_filter:
                continue
            bloom_filter.add(kmer)
            # unique_kmers_set.add(kmer)
            kmers_added += 1

        total_seq += 1

    return total_kmers, bloom_filter, total_seq, kmers_added

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
    parser.add_argument('-n', '--cpus', type=int, default=0.01,
                    help='The number of cores to use (default 1)')

    args = parser.parse_args()
    
    '''if choose to build bloom filter, 
       load exclude bloom filter, 
       choose appropriate kmer length (dafault 31) and false positive, 
       build new bloom filter
    '''
    
    if args.data_structure=='bloom':
        if args.exclude_filter:
            print("Loading exclude bloom filter...")
            exclude_filter_bf = BloomFilter.load(args.exclude_filter)
            k = exclude_filter_bf.kmer_length
            print(f"Using k-mer length {k} from the exclude bloom filter.")
        else:
            exclude_filter_bf = None
            k = args.kmer_length

        if args.cpus > 1:
            print("Since using parallel processing, loading the whole file into memory...")
            records = list(SeqIO.parse(args.contamination_fasta, "fasta"))
            print("finished loading sequences")

        if args.expected_elements:
            n_unique = args.expected_elements
        else:
            #estimate number of unique kmer and total kmer

            if args.cpus > 1:
                n_unique = estimate_unique_kmers_parallel(records, k, args.cpus, exclude_filter_bf if args.exclude_filter else None)
            else:
                n_unique = estimate_unique_kmers(args.contamination_fasta, k, exclude_filter_bf if args.exclude_filter else None)
            
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

        if args.cpus > 1:
            print("Building Bloom filter with parallel processing...")
            total_kmers = 0
            total_seq = 0
            not_excluded = 0

            # Distribute records for the workers in a round-robin fashion
            chunks = [[] for _ in range(args.cpus)]
            for i, record in enumerate(records):
                chunks[i % args.cpus].append(record)

            with Pool(processes=24, initializer=init_bf_worker, initargs=(exclude_filter_bf, n_unique, false_positive_rate, k)) as pool:
                with tqdm(total=len(chunks), desc="Processes adding k-mers") as pbar:
                    func = partial(process_sequence_chunk, k=k)
                    for result in pool.imap_unordered(func, chunks):
                        total_kmers += result[0]
                        bloom_filter.combine(result[1])
                        total_seq += result[2]
                        not_excluded += result[3]

                        pbar.update(1)

            if total_kmers > 0:
                print(f"Total kmer not excluded: {not_excluded}")
                print(f"Total k-mers processed: {total_kmers}")
                print(f"Total seq processed: {total_seq}")
            else:
                print("No k-mers were processed.")
        else:
            print("Building Bloom filter...")
            total_kmers = 0
            kmers_added = 0
            #use SeqIO get each sequence from fasta file
            for record in tqdm(SeqIO.parse(args.contamination_fasta, "fasta"), desc="Building Bloom filter"):
                seq = str(record.seq).upper()
                for kmer in generate_kmers(seq, k):
                    total_kmers += 1
                    if args.exclude_filter and kmer in exclude_filter_bf:
                        continue
                    bloom_filter.add(kmer)#add kmer to bloom filter
                    kmers_added += 1 #count the number of unique kmers
            if total_kmers > 0:
                print(f"{kmers_added} out of {total_kmers} total kmers were added.")
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
            exclude_filter = BCuckooFilter.load(args.exclude_filter)
            k = exclude_filter.kmer_length
            print(f"Using k-mer length {k} from the exclude CuckooFilter.")
        else:
            k = args.kmer_length
            print(f"Using k-mer length of {k}")

        if args.expected_elements:
            n_unique = args.expected_elements
        else:
            #estimate the number of unique kmer and total kmer
            n_unique = estimate_unique_kmers(args.contamination_fasta, k, exclude_filter if args.exclude_filter else None)

        '''calculate the bucket size and estimate the size of fingerprint and capacity
            the equation get from https://stackoverflow.com/questions/57555236/how-to-size-a-cuckoo-filter
        '''
        if args.false_positive_rate<0.002:
            bucket_size=4
        else:
            bucket_size=2
        
        print(f"Using bucket size of {bucket_size}")

        if args.capacity_of_cuckoofilter:
            capacity=args.capacity_of_cuckoofilter
        else:
            # fingerprint_size=math.ceil((math.log2(1 / args.false_positive_rate)+math .log2(2*bucket_size)))
            # print(f"Using fingerprint size of {fingerprint_size}")
            if bucket_size==4:
                capacity=math.ceil((n_unique/0.95/bucket_size))
            else:
                capacity=math.ceil((n_unique/0.84/bucket_size))
            print(f"Using capacity size of {capacity}")
        #build new cuckoo filter
        cuckoo = BCuckooFilter(capacity,args.false_positive_rate,k, bucket_size)
        # cuckoo_size_bytes = cuckoo.__sizeof__()
        print(f"size: {len(cuckoo.buckets) / 8 / (1024 ** 3):.4f} GB")
        print("Building Cuckoo Filter...")

        total_kmers = 0
        kmers_added = 0
        #use SeqIO get the sequence from fasta file
        for record in tqdm(SeqIO.parse(args.contamination_fasta, "fasta"), desc="Building Cuckoo filter"):
            seq = str(record.seq).upper()
            for kmer in generate_kmers(seq, k):
                total_kmers += 1#count the number of total kmer
                if args.exclude_filter and kmer in exclude_filter:
                    continue
                cuckoo.insert(kmer)#add kmer in cuckoo filter
                kmers_added += 1#count the number of unique kmer
        if total_kmers > 0:
            print(f"{kmers_added} out of {total_kmers} total kmers were added.")
        else:
            print("No k-mers were processed.")
        cuckoo.save(args.output_filter)#save cuckoo filter as output file
        print(f"cuckoo filter saved to {args.output_filter}")

if __name__ == "__main__":
    main()
