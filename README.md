[![PyPI version](https://img.shields.io/pypi/v/KmerDecon.svg)](https://pypi.org/project/KmerDecon/)
# KmerDecon

KmerDecon is a fast, memory-efficient tool for decontaminating sequencing reads using Bloom filters or Cuckoo filters. It generate detailed reports of contaminants in sequencing data.

## Authors
- Yujia Feng
- Xiaoyi Chen
- Yuxiang Li

## Installation

### Prerequisites:

- Python 3.6 or higher
- pip package manager

### Install:

Run the following command:
  ```
  pip install KmerDecon
  ```

## Usage

### 1. Building the Bloom Filter or Cuckoo Filter

Generate a Bloom filter from contamination source sequences. Generate a Cuckoo filter, use `-s cuckoo`. Use `kbuild --help` for more detail.

```
kbuild -c contamination.fasta -s bloom -o contamination_filter.bf
```


**Optional Arguments:**

- `kmer-length`: Length of k-mers to generate (e.g., 31). If not provided, the tool determines the optimal k-mer length automatically.
- `expected-elements`: Expected number of unique k-mers. If not provided, it is estimated using HyperLogLog.
- `exclude-filter`: A .bf filter or .cms file path. If provided, any k-mers present in the excluded filter will not be encoded into the new build filter.
- `max-memory`: Maximum memory in GB for the Bloom filter. Adjusts parameters to fit within this limit.
- `false-positive-rate`: Desired false positive rate (default: 0.01).

if choose build Cuckoo filter:
- `capacity-of-cuckoofilter`: The capacity of cuckoo filter

### 2. Decontaminating Reads

Filter out contaminated reads from your sequencing data. Use `kdecon --help` for more detail.

Use bloom filter:
```
kdecon -i reads.fastq -d example_filter/hg38.bf -s bloom -o output
```
Use `-s cuckoo` for Cuckoo filter.

**Optional Arguments:**

- `threshold`: Fraction of matching k-mers to consider a read contaminated (default: 0.5).
- `kmer-length`: Length of k-mers used. If not provided, the k-mer length from the Bloom filter is used.
- `mode`: Operation mode, either filter (default) or states.
  - filter: Filters reads based on contamination levels.
  - states: Generates a states.csv report with contamination statistics. Columns:
	- {filter}_avgSimilarity: The average fraction of matching k-mers across all reads in that file for each filter.
	- {filter}_percentReadsPassing: The percentage of reads passing the threshold for each filter.

## Performance

### Highlights

- With default parameters, we achieves FPR = 0.002%, FNR = 0.05% on simulated human reads decontamination task.
- KmerDecon is memory efficient and uses 10 bits / kmer. (Popular too Kraken2 uses 32 bits / kmer)
- KmerDecon is fast and takse 5 min to filter 1 million reads of 150bp each (kraken2 takes ~8min, both on single thread)
- Multi-threads parallel building supported.

### Full Reports
- To read the full performance report, please see: [Here](https://drive.google.com/file/d/1shEp8LZAC5w_qR0p8BzsjBeAsOffvDIH/view?usp=sharing)
- To recreate the results on the report, please see: [Here](reproducibility.md)


## Dependencies

- `bitarray>=2.1.0`
- `biopython>=1.78`
- `mmh3>=2.5.1`
- `hyperloglog>=0.0.12`

Install dependencies manually with:

```bash
pip install -r requirements.txt
```
## Referenced Code
 The python module of cuckoofilter is adapted from:

 Author: Huy Do

 Repository: https://github.com/huydhn/cuckoo-filter/blob/master/cuckoo/filter.py
 
 License: MIT

## Contributing

Contributions and PRs are welcome!

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For questions or suggestions, please open an issue.
