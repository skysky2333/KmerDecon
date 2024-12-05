[![PyPI version](https://img.shields.io/pypi/v/KmerDecon.svg)](https://pypi.org/project/KmerDecon/)
# KmerDecon

KmerDecon is a fast, memory-efficient tool for decontaminating sequencing reads using Bloom filters or Cuckoo filters. It generate detailed reports of contaminants in sequencing data.

## Authors
- Yujia Feng
- Xiaoyi Chen
- Yuxiang Li


## Features

- **Automatic Parameter Optimization**: Automatically determines the optimal k-mer length and adjusts parameters based on desired memory and false positive rate using tools like HyperLogLog.
- **Speed**: Utilizes efficient hashing with MurmurHash3 for fast k-mer processing.
- **Memory Efficiency**: Employs Bloom filters with dynamic sizing to balance memory usage and accuracy, capable of handling billions of k-mers with minimal RAM.
- **Scalability**: Suitable for large datasets, such as whole-genome sequencing reads and large contamination sources like the human genome.
- **Detailed Reporting**: Generates comprehensive reports on contamination levels across multiple samples and filters.
- **Real-Time Processing**: Allows for decontamination during data streaming or generation, providing immediate feedback and contaminant removal. (TODO)

## Installation

### Prerequisites:

- Python 3.6 or higher
- pip package manager

### Steps:

1. Install directory:
   ```
   pip install KmerDecon
   ```

2. Alternatively, to get the lastest version, you can clone the repository:

    ```
    git clone https://github.com/skysky2333/KmerDecon
    cd KmerDecon
    pip install .
    ```

## Usage

### 1. Building the Bloom Filter or CountMin sketches

Generate a Bloom filter from contamination source sequences. Use `kbuild --help` for more detail.

```
kbuild -c contamination.fasta -s bloom -o contamination_filter.bf
```
Generate a CountMin sketches from contamination source sequences. Use `kbuild --help` for more detail.

```
kbuild -c contamination.fasta -s cuckoo -o contamination_filter.cf 
```


**Optional Arguments:**

- `kmer-length`: Length of k-mers to generate (e.g., 31). If not provided, the tool determines the optimal k-mer length automatically.
- `expected-elements`: Expected number of unique k-mers. If not provided, it is estimated using HyperLogLog.
- `exclude-filter`: A .bf filter or .cms file path. If provided, any k-mers present in the excluded filter will not be encoded into the new build filter.
- `max-memory`: Maximum memory in GB for the Bloom filter. Adjusts parameters to fit within this limit.
- `false-positive-rate`: Desired false positive rate (default: 0.001).

if choose build Cuckoo filter:
- `capacity-of-cuckoofilter`: The capacity of cuckoo filter
- `fingerprint-size-of-cuckoofilter`: The depth of cuckoo filter

### 2. Decontaminating Reads

Filter out contaminated reads from your sequencing data. Use `kdecon --help` for more detail.

Use bloom filter:
```
kdecon -i reads.fastq -d example_filter/hg38.bf -s bloom -o output
```
Use countmin sketch:
```
kdecon -i reads.fastq -d example_filter/hg38.cf -s cuckoo -o output
```

**Optional Arguments:**

- `threshold`: Fraction of matching k-mers to consider a read contaminated (default: 0.5).
- `kmer-length`: Length of k-mers used. If not provided, the k-mer length from the Bloom filter is used.
- `mode`: Operation mode, either filter (default) or states.
  - filter: Filters reads based on contamination levels.
  - states: Generates a states.csv report with contamination statistics. Columns:
	- {filter}_avgSimilarity: The average fraction of matching k-mers across all reads in that file for each filter.
	- {filter}_percentReadsPassing: The percentage of reads passing the threshold for each filter.


## Dependencies

- `bitarray>=2.1.0`
- `biopython>=1.78`
- `mmh3>=2.5.1`
- `hyperloglog>=0.0.12`

Install dependencies with:

```bash
pip install -r requirements.txt
```
## Referenced Code
 The python module of cuckoofilter and bucketis are adapted from:
 Author: Michael The
 Repository: https://github.com/michael-the1/python-cuckoo/tree/master
 License: MIT

## Contributing

Contributions and PRs are welcome!

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For questions or suggestions, please open an issue.
