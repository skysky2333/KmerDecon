# KmerDecon

KmerDecon is a fast, memory-efficient tool for decontaminating sequencing reads using Bloom filters. It allows for on-the-fly filtering of contaminants in sequencing data.
## Authors
Yujia Feng, Xiaoyi Chen, Yuxiang Li

## Features

- **Speed**: Utilizes Bloom filters for constant-time membership queries.
- **Memory Efficiency**: Compact representation reduces RAM usage.
- **Real-Time Processing**: Decontamination during data streaming or generation.

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
build-bloom-filter --contamination-fasta contamination.fasta --kmer-length 31 --output-filter contamination_filter.bf
```

**Arguments:**

- `--contamination-fasta`: Path to the FASTA file containing contamination sequences.
- `--kmer-length`: Length of k-mers to generate (e.g., 31).
- `--output-filter`: Output filename for the Bloom filter.

### 2. Decontaminating Reads

Filter out contaminated reads from your sequencing data.

```bash
decontaminate-reads --input-reads reads.fastq --bloom-filter contamination_filter.bf --threshold 0.5 --kmer-length 31 --output-reads decontaminated_reads.fastq
```

**Arguments:**

- `--input-reads`: Path to the input FASTQ file with sequencing reads.
- `--bloom-filter`: Path to the Bloom filter file generated earlier.
- `--threshold`: Fraction of matching k-mers to consider a read contaminated (default: 0.5).
- `--kmer-length`: Length of k-mers used (must match the value used when building the filter).
- `--output-reads`: Output FASTQ file for decontaminated reads.

### Examples

**Building the Bloom Filter**

```bash
build-bloom-filter --contamination-fasta human_genome.fasta --kmer-length 31 --output-filter human_bloom_filter.bf
```

**Decontaminating Reads**

```bash
decontaminate-reads --input-reads sample_reads.fastq --bloom-filter human_bloom_filter.bf --threshold 0.5 --kmer-length 31 --output-reads decontaminated_reads.fastq
```

## Dependencies

- `bitarray>=2.1.0`
- `biopython>=1.78`

Install dependencies with:

```bash
pip install -r requirements.txt
```

## Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Commit your changes.
4. Submit a pull request with a detailed description.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

For questions or suggestions, please open an issue.
