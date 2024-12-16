# Recreating the Results on the Report

All the instructions were done on a Linux computer.

## 1. Download NCBI datasets command line tool
The NCBI datasets command line tool can be downloaded [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/), or alternatively run the follow commands.
```
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
chmod +x datasets
```

## 2. Downloading the Required Genomes

We created 5 fasta files. One for the human genome, one for the mouse genome, one for the genomes for the genuses of (Arthrobacter, Burkholderia, Chryseobacterium, Ochrobactrum, Pseudomonas, Ralstonia, Rhodococcus, Sphingomonas, Corynebacterium, Propionibacterium and Streptococcus) which we refer to as common contaminants, one for all bacteria, and one for fruit fly.

**Download human genome:**
```
./datasets download genome taxon human --reference --filename human.zip
```

**Download mouse genome:**
```
./datasets download genome taxon "mus musculus" --reference --filename mouse.zip
```

**Download common contaminants genomes:**  
You can run the command line for each of the genus similar to how we did it for the human and mouse genome, but with an extra flag of `--assembly-level chromosome,complete`. Alternatively, you can use the `download-common-contam.py` script to download the genomes instead of manual typing each command.  

Then we simply combined each of the fasta files into one. You can also do this, by using the `combine.py` script, after moving all the folders into one directory. Remember to change the `base_dir` and `output_file` variable to the folder/filename you are using/want.

**Download all bacteria genomes:**  
WARNIG: The combined size exceeds 20GB. These genomes were only used to generate simulated reads and was not used to build any data structures. If you do not wish to download this, we can provide you with the reads that we generated.

```
./datasets download genome taxon 2 --reference --assembly-level chromosome,complete --filename bacteria.zip
```

The downloaded genomes are then combined into one fasta file. The same `combine.py` script can be used to combine them.

**Download fruit fly genome:**
```
./datasets download genome taxon 7227 --reference --assembly-level chromosome,complete --filename fruitfly.zip
```


## 3. Building the Bloom filters

WARNING: While building the Bloom filters using parallel processing significantly reduces the time, it requires significantly more memory as each worker maintains its own Bloom filter and sequences. For reference, constructing the human Bloom filter using 24 processes has around 150GB peak usage of memory.

**TIME REQUIRED TO BUILD**

The following time was obtained when running the program on the gradx computer at 24 cores.

For the human genome:
- 30min for HLL
- 55min to build filter

For the common contaminants:
- 25min for HLL
- 35min to build filter

For the mouse genome:
- 40min for HLL
- 60min to build filter

For the fruit fly genome:
- 4min for HLL
- 7min to build filter

The HLL step can be skipped by providing the algorithm with the expected number of elements with the flag `-e`. We will provide this number for each of the Filters since we have already ran the program and know the number. Use to `-n` to indicate the number of cores to use if you want to build it using parallel processing.


**First build the human Bloom filter**

```
kbuild -s bloom -c <human genome> -k 31 -o <output name>
```

To skip the HLL step, include the flag `-e 2562928088`.

**Build the common contaminants Bloom filter excluding Human**

```
kbuild -s bloom -c <commmon contaminants genome> -x <human Bloom filter> -k 31 -o <output name>
```

The `<human Bloom filter>` should be the Bloom filter file generated from the previous step. To skip the HLL step, include the flag `-e 1809667235`.

**Build the mouse Bloom filter excluding Human**

```
kbuild -s bloom -c <mouse genome> -x <human Bloom filter> -k 31 -o <output name>
```

To skip the HLL step, include the flag `-e 2203333902`.

**Build the fruit fly Bloom filter**
```
kbuild -s bloom -c <fruit fly genome> -k 31 -o <output name>
```

To skip the HLL step, include the flag `-e 126346179`.

## 4. Build the Cuckoo filter

The cuckoo filter takes significantly longer to build compared to Bloom filters. The build time for the fruit fly Cuckoo filter took around 30 minutes. To build it run the following code:

```
kbuild -s cuckoo -c <fruit fly genome> -k 31 -o <output name>
```

## 5. Generating simulated reads

The simulated reads were generated using [InSilicoSeq](https://insilicoseq.readthedocs.io/en/latest/#). For each of the fasta files (5 in total), we generated 2 million reads only half of which were used (the files ending with R1) as it generates paired-end reads.

**Install the tool**
```
pip install InSilicoSeq
```

**Generating the reads**
```
iss generate --genomes <fasta file> --model novaseq --n_reads 2M --output <output name>
```

After generating the reads, two files will be created, ...R1.fastq and ...R2.fastq, only use the R1 file for our program.

## 6. Decontaminating reads

Use the `--mode states` to obtain the statistic described in our report.

```
kdecon -i <fastq reads> -d <Bloom filter> -s bloom -o <output dir> --mode states
```

**Decontamination that we preformed**

human Bloom filter:
- human reads
- mouse reads
- common contaminants reads
- bacteria reads

common contaminants Bloom filter:
- human reads
- common contaminants reads
- bacteria reads

mouse Bloom filter:
- human reads
- mouse reads

fly Bloom filter:
- fly reads

fly Cuckoo filter:
- fly reads