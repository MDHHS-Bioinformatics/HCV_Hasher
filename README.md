## Introduction

**HCV-Hasher** is a bioinformatics python script designed for the analysis of WGS Hepatitis C haplotypes to determine possible linkages between samples. The script is built using Python:3.10.16 and the following modules are required for the script to function :
- `mmh3==5.1.0` 
- `pandas==2.2.3` 
- `biopython==1.85` 
- `numpy==2.2.4`

HCV-Hasher was inspired by principals from [CDC's GHOST](https://pmc.ncbi.nlm.nih.gov/articles/PMC5731493/) and [Kraken2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) to be able to identify possible linkages between genetically similar samples. The algorithm generates k-mers for each unique haplotype in a sample. A hash table is used to efficiently store and retrieve associated k-mer/sample haplotype information. The Jacquard index is the equation used to calculate the percent similarity for all haplotypes between all pairs of samples. 

## Input/Output Options

```bash
usage: hcv_hasher.py [-h] [--table_size TABLE_SIZE] --mode {initialize,update,load} [--samples_dir_to_add SAMPLES_DIR_TO_ADD] [--samples_dir_to_compare SAMPLES_DIR_TO_COMPARE] [--previous_hash_table PREVIOUS_HASH_TABLE] [--new_hash_table_name NEW_HASH_TABLE_NAME] [--kmer_size KMER_SIZE] [--save_table] [--linkage_threshold LINKAGE_THRESHOLD][--version]

Hash table and kmer-ization implementation for HCV clustering

options:
  -h, --help            show this help message and exit
  --table_size TABLE_SIZE
                        Size of the hash table to initialize with
  --mode {initialize,update,load}
                        Mode to perform: initialize, update, or load
  --samples_dir_to_add SAMPLES_DIR_TO_ADD
                        Directory of samples to add to the hash table
  --samples_dir_to_compare SAMPLES_DIR_TO_COMPARE
                        Directory of samples to compare to the hash table
  --previous_hash_table PREVIOUS_HASH_TABLE
                        Path to a previous hash table to load in and use
  --new_hash_table_name NEW_HASH_TABLE_NAME
                        Name of the final outputted Hash table. If a hash table is being updated, the updated table will be saved to the new name
  --kmer_size KMER_SIZE
                        Size of the k-mer to be used in the analysis
  --save_table          Save the hash table after creation
  --linkage_threshold LINKAGE_THRESHOLD
                        Percent similarity threshold to store the specific haplotype linkages of potentially linked samples. Default:50
  --version             show program's version number and exit
```

## Modes

HCV Hasher has three main options for `--mode`: 
- `initialize` 
- `update` 
- `load` 
#### `initialize`
The  `initialize` option creates a brand new hash table. The only true required option for this mode is `--samples_dir_to_add`, in order to supply the path to a directory of fasta files, where each fasta file is one sample containing all its unique haplotypes. 
#### `update`
The `update` option is designed to be used when a previously created hash table needs to be updated with samples that are *not* in the current hash table. Once again, `--samples_dir_to_add` is needed with the directory path to add new samples to the hash table. Make sure that the `--kmer_size` used in your previous hash table is the same across updates. If you would also like to compare these new samples at the same time, you can give the same directory path to `--samples_dir_to_compare` in your command.
#### `load`
The `load` option is used to load a previously created hash table, with `--previous_hash_table`. This is primarily used to compare samples, with `--samples_dir_to_compare` providing a directory to a set of fasta files.

### Outputs
When using the `--samples_dir_to_compare` option, you can expect the generation of the following files
- `within_sample_percent_similarities.csv`
- `between_sample_percent_similarities.csv`
- `haplotype_linkages.csv`

For `within_sample_percent_similarities.csv` you can expect the following output format: 
| sample_1 | sample_2 | percent_similarity |
|--------|---------|---------|
| SAMPLE-1_0 | SAMPLE-1_1 | 97.8% |
| SAMPLE-1_0 | SAMPLE-1_2 | 97.6% |
| SAMPLE-1_0 | SAMPLE-1_3 | 97.1% |
| SAMPLE-1_1 | SAMPLE-1_2 | 96.5% |
Every row in the file contains the comparison of two haplotypes from the *same* sample and their percent similarities. Every possible haplotype pair comparison per sample is included in this file 

For `between_sample_percent_similarities.csv` you can expect the following output format: 
| sample_1 | sample_2 | percent_similarity |
|--------|---------|---------|
| SAMPLE-1 | SAMPLE-2 | 82.8% |
| SAMPLE-1 | SAMPLE-3 | 67.6% |
| SAMPLE-4 | SAMPLE-5 | 30.1% |
| SAMPLE-6 | SAMPLE-8 | 27.5% |
Every row in the file contains the *highest* percent similarity result after comparing all possible haplotype pairs between a pair of samples. Every fasta present in the directory specified by `--samples_dir_to_compare` will have a result for every sample present in the hash table. Comparisons of the same sample are not performed. 

For `haplotype_linkages.csv` you can expect the following output format: 
| sample_1 | sample_2 | percent_similarity |
|--------|---------|---------|
| SAMPLE-1_4 | SAMPLE-2_1 | 82.8% |
| SAMPLE-1_11 | SAMPLE-3_9 | 67.6% |
This file will contain similar data to `between_sample_percent_similarities.csv`, the difference however is that instead of just the sample name being present for each row, the sample name alongside the specific haplotype will be in each row. The format will be `SAMPLE-NAME_HAPLOTYPE#` Additionally, this file is filtered down to a desired threshold using `--linkage_threshold` to only provide the most genetically similar samples. 

### Input Fasta Formatting
When adding or comparing sequence data (with `--samples_dir_to_add` or `--samples_dir_to_compare`) you must specify a directory containing your fasta file(s). Fasta files names should be the name of the sample followed by the genotype of the sample using a '-'. So for a list of samples: `EX-12345`, `EX-23456`, and `EX-34567`, we would append their genotypes at the end to look like : `EX-12345-1a`, `EX-23456-2b`, and `EX-34567-3a`. The contents of each individual fasta file should be in a multi-fasta format where each sequence is a unique haplotype for that sample. The format should be `>SAMPLE-NAME-GENOTYPE_HAPLOTYPE#`: 
```bash
>EX-12345-1a_0
ACCCTTGGCAGGCA....
>EX-12345-1a_1
ACCCTTGGCAGGCG....
>EX-12345-1a_2
ACCCTTGGCAGGCC....
>EX-12345-1a_3
ACCCTTGGCAGGCT....
```
The genotype should always be specified in each fasta as multi-genotype samples can then be analyzed separately per unique genotype. 

## Script Use Examples
You first need to create your initial database. You would start with running the following command and adjusting any parameters you see as necessary

```bash
  python hcv_hasher.py \
  --mode initialize \
  --samples_dir_to_add /path/to/fastas/ \
  --table_size 100000 \
  --save_table \
  --new_hash_table_name hcv_hash_table
```

Now that you have a created table, you may need to update it with new sequences, to do so, run a command like the following:
```bash
  python hcv_hasher.py \
  --mode update \
  --samples_dir_to_add /path/to/extra/fastas/ \
  --previous_hash_table hcv_hash_table.pkl
  --save_table \
  --new_hash_table_name hcv_hash_table
```
Now that you have a hash table with your sequences, you likely want to be able to see how similar your samples are. To do so, run a command like the following: 
```bash
  python hcv_hasher.py \
  --mode load \
  --samples_dir_to_compare /path/to/extra/fastas/ \
  --previous_hash_table hcv_hash_table.pkl
  --linkage_threshold 60
```


### A Note on Percent Similarity Interpretation and Relation to K-mer size Selection

Traditional taxonomic identification tools such as Kraken2 use K-mers that are large enough to ensure species-level specificity but not large enough to where sensitivity is lost. Previous development versions of this script used a default K-mer size of 25, which is closer to the default K-mer size of 35 used by Kraken2. Intuitively, since our analysis are purely within a single species, we are looking for smaller levels of differentiation, meaning that a smaller K-mer size could be used. Based on initial performance, k=25 was satisfactory for being able to create distinct distributions to examine within-sample, same-genotype, and different-genotype. 

Mathematically, there [*technically*](https://pmc.ncbi.nlm.nih.gov/articles/PMC11874746/#GR279452JENC11:~:text=and%20genome%20profiling.-,Essential%20properties%20of%20k%2Dmers,-The%20choice%20of) is an optimal K-mer size that could be employed so that most k-mers in the genome will be found only once (i.e. having only unique k-mers and no repeats). To estimate the minium k-mer size for this, we can compute the expected number of occurences a given k-mer occurs in genome of length $G$ (HCV genome is ~9.6kb). 
- This estimation is:
$G / 4^k$
  - Where 4 is the alphabet size (i.e. the possible expected base of A,C,G,T), G is the genome length in bps, and k is the optimial k-mer size we are looking for
- Which can be simplified to:
$\log_4(G) \geq k$
- For the HCV genome this calculation would be:
$9600 / 4^k$
- $\log_4(9600) \geq k$
- $7 \geq k$

Therefore, the theoretical minimum k-mer size of 7 should be chosen for HCV to ensure that each k-mer is unique for HCV.Practically however, choosing a k-mer size that is slightly larger than the theoretical minimum can be advantageous to further mitigate the effects of repetitive regions, gaps ('---'), as well as incomplete genomes (i.e. long stretches of Ns).

Further evaluation of the k-mer size was performed by dropping the k-mer size in increments of five: 25, 20, 15, 10, and the theoretical minimum size of 7. Adjusting k-mer sizes in this range changes the distributions of percent similarities when looking across the three main distributions of  within-sample, same-genotype, and different-genotype. Dropping the k-mer size in this range shifts the distributions upwards (higher percent similarity). The possible linkages don't actually change as distributions shift. The thresholds for what should be considered a linkage must be adjusted as the k-mer size changes.

For now, our testing leans towards using a k-mer size of 10 for optimal linkage analysis. At k=10, the percentage of unique 10-mers that occur only once (per sample across all haplotypes) is >99%, and that less than <1% of 10-mers occur more than once. 