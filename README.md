# KmerDecon

KmerDecon is a fast, memory-efficient tool for decontaminating sequencing reads using Bloom filters. It generate detailed reports of contaminants in sequencing data.
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

1. Clone the repository:

    ```bash
    git clone https://github.com/skysky2333/KmerDecon
    ```

2. Navigate to the project directory:

    ```bash
    cd KmerDecon
    ```

3. Install the package:

    ```bash
    pip install .
    ```

## Usage

### 1. Building the Bloom Filter

Generate a Bloom filter from contamination source sequences.

```bash
build-bloom-filter --contamination-fasta contamination.fasta --output-filter contamination_filter.bf
```

**Optional Arguments:**

- `kmer-length`: Length of k-mers to generate (e.g., 31). If not provided, the tool determines the optimal k-mer length automatically.
- `max-memory`: Maximum memory in GB for the Bloom filter. Adjusts parameters to fit within this limit.
- `false-positive-rate`: Desired false positive rate (default: 0.001).
- `expected-elements`: Expected number of unique k-mers. If not provided, it is estimated using HyperLogLog.

### 2. Decontaminating Reads

Filter out contaminated reads from your sequencing data.

```bash
decontaminate-reads --input-reads reads.fastq(can_also_be_directory) --bloom-filter contamination_filter.bf(can_also_be_directory) --output-dir output_directory
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

## Contributing

Contributions and PRs are welcome!

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For questions or suggestions, please open an issue.
