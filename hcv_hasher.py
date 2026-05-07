#!/usr/bin/env python
# -*- coding: utf-8 -*-
# hcv_hasher.py
# A script to implement a hash table for storing k-mer information sample IDs

__version__ = "1.0.3"

import pandas as pd 
from Bio import AlignIO
from Bio import SeqIO
import os 
import pickle
import argparse
import logging
from collections import Counter
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

 

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class FullHashTable:
    """
    Dictionary-backed k-mer index.
    Stores:
        {encoded_kmer_int : set(sample_ids)}
    """

    def __init__(self):
        self.table = defaultdict(set)

    def insert(self, key, sample_id):
        """
        Insert a k-mer and associated sample ID.
        """
        self.table[key].add(sample_id)

    def search(self, key):
        """
        Return set of sample IDs associated with k-mer.
        """
        return self.table.get(key, None)

    def count_filled(self):
        """
        Count number of unique kmers stored.
        """
        filled = len(self.table)
        empty = 0
        return filled, empty

    def save_table(self, filename):
        """
        Save hash table to file.
        """
        with open(filename, "wb") as f:
            pickle.dump(dict(self.table), f)
        print(f"Hash table saved to {filename}")

    @staticmethod
    def load_table(filename):
        """
        Load hash table from file.
        Only supports dict-based format.
        """
        with open(filename, "rb") as f:
            table = pickle.load(f)

        if not isinstance(table, dict):
            raise ValueError(
                "Unsupported table format. Expected dict-based table."
            )

        new_ht = FullHashTable()
        new_ht.table = defaultdict(set, table)

        print(f"Hash table loaded from {filename}")
        return new_ht

    def get_sample_id_counts(self):
        """
        Count number of unique kmers per sample ID.
        """
        sample_id_counts = {}

        for sample_ids in self.table.values():
            for sample_id in sample_ids:
                sample_id_counts[sample_id] = (
                    sample_id_counts.get(sample_id, 0) + 1
                )

        return sample_id_counts

    def get_sample_id_counts_df(self):
        """
        Return DataFrame with sample ID counts.
        """
        sample_id_counts = self.get_sample_id_counts()
        return pd.DataFrame(
            list(sample_id_counts.items()),
            columns=["sample_id", "counts"]
        )

    def display_table(self):
        """
        Print table contents.
        """
        print(dict(self.table))


def iterate_over_samples(samples_dir, hash_table,kmer_size):
    """
    Iterates over all sample files in the given directory, processes each file to extract k-mers,
    and inserts them into the provided hash table. Resizes the hash table if the load factor exceeds 0.65.

    Args:
        samples_dir (str): The directory containing sample files in FASTA format.
        hash_table (HashTable): The hash table to insert k-mers into.

    Returns:
        HashTable: The updated hash table with k-mers from all processed samples.

    Raises:
        OSError: If there is an issue reading files from the samples directory.
        ValueError: If there is an issue with the format of the sample files.
    """
    for filename in os.listdir(samples_dir):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            filepath = os.path.join(samples_dir, filename)
            records = list(SeqIO.parse(filepath, "fasta"))
            print(f"Processing {filename} with {len(records)} sequences.")
            for record in records:
                sample_name = record.id
                if not sample_name.rsplit('_', 1)[-1].isdigit():
                    raise ValueError(f"FASTA header '{sample_name}' in '{filename}' does not end with '_N', where N is the haplotype number")
                #track seen kmers to prevent double counting
                seen_kmers = set()
                #add each unique kmer to the hash table
                for kmer in iter_kmers(str(record.seq), kmer_size):
                    if kmer in seen_kmers:
                        continue
                    seen_kmers.add(kmer)
                    hash_table.insert(kmer, sample_name)
    print("\nAll samples processed and added to hash table")
    return hash_table

def iter_kmers(seq, k):
    """
    Rolling 3-bit integer encoder.
    Keeps ambiguous bases as part of k-mers.
    No resets.
    """
    BASE3BIT = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "N": 4,
    "0": 4,
    "-": 5
    }
    mask = (1 << (3 * k)) - 1
    val = 0
    length = 0

    for base in seq.upper():
        code = BASE3BIT.get(base, 4)  # unknowns → N

        val = ((val << 3) | code) & mask
        length += 1

        if length >= k:
            yield val


GLOBAL_HASH_TABLE = None
GLOBAL_ID_COUNTS = None

def init_worker(hash_table, id_counts):
    global GLOBAL_HASH_TABLE, GLOBAL_ID_COUNTS
    GLOBAL_HASH_TABLE = hash_table
    GLOBAL_ID_COUNTS = id_counts
def process_single_fasta(args):
    """
    Iterate overa single fasta file and compute the within sample and between sample similiarties for each sequence/haplotype in the fasta

    Args:
        (args): A tuple containing the following elements:
            fasta (str): The name of the fasta file to process.
            samples_dir (str): The directory containing the fasta file.
            hash_table (FullHashTable): The hash table to search for k-mers.
            kmer_size (int): The size of the k-mers to be used in the analysis.
            linkage_threshold (int): The percent similarity threshold to store specific haplotype linkages.
            id_counts (dict): A dictionary containing the count of unique k-mers for each sample stored in the FullHashTable.

    Returns:
        tuple: A tuple containing three dictionaries:
            - within_sample_dict (dict): A dictionary with keys as tuples of sample pairs and values as their percent similarity for within-sample comparisons.
            - between_sample_dict (dict): A dictionary with keys as tuples of sample pairs and values as their percent similarity for between-sample comparisons.
            - haplotype_linkage_dict (dict): A dictionary with keys as tuples of sample pairs and values as tuples containing the best haplotype pair and their percent similarity for potentially linked samples.
    """
    # unpack the arguments passed from the the multiprocessing call
    (fasta, samples_dir, kmer_size, linkage_threshold) = args
    

    # initialize results dictionaries
    within_sample_dict = {}
    between_sample_dict = {}
    haplotype_linkage_dict = {}
    #format sample name and path properly
    fasta_path = os.path.join(samples_dir, fasta)
    #iterate through each haplotype/sequence in the fasta file
    for record in SeqIO.parse(fasta_path, "fasta"):
        sample_name = record.id
        if not sample_name.rsplit('_', 1)[-1].isdigit():
            raise ValueError(f"FASTA header '{sample_name}' in '{fasta}' does not end with '_N', where N is the haplotype number")
        primary_base = sample_name.rsplit('_', 1)[0]
        #store the shared number of kmers with other haplotypes
        temp_counter = Counter()
        #track seen kmers to prevent double counting
        seen_kmers = set()
        #start unique kmers counts before iterating
        kmer_count = 0
        #iterate through the kmers
        for kmer in iter_kmers(str(record.seq), kmer_size):
            #if the kmer has already been counted/seen, skip it
            if kmer in seen_kmers:
                continue
            seen_kmers.add(kmer)
            kmer_count += 1
            #find the sample that share the current kmer
            result_set = GLOBAL_HASH_TABLE.search(kmer)
            #for each sample, count shared kmers
            if result_set:
                temp_counter.update(result_set)
        #calculate similarities for each matched haplotype
        for hit_sample, hit_counts in temp_counter.items():
            #get the base sample name for hit sample
            hit_base = hit_sample.rsplit('_',1)[0]
            #get the number of unique kmers in the matched haplotype
            hit_sample_kmer_count = GLOBAL_ID_COUNTS.get(hit_sample,0)
            #calculate the union size of unique kmers between the two haplotypes
            total_unique_counts = (kmer_count + hit_sample_kmer_count - hit_counts)
            #calculate the percent similiarity (jaccard index)
            percent_similarity = round((hit_counts / total_unique_counts) * 100, 3)
            ### WITHIN SAMPLE ###
            if primary_base == hit_base:
                #skip self comparisons
                if sample_name == hit_sample:
                    continue
                #Create consistent ordering of pairs
                key = tuple(sorted([sample_name, hit_sample]))
                within_sample_dict[key] = percent_similarity
            ### BETWEEN SAMPLE ###
            else:
                #create sample pair alphabetically
                sample_pair_key = (primary_base,hit_base
                ) if primary_base < hit_base else (hit_base,primary_base)
                #only keep the highest percent similiarity between two samples
                if percent_similarity > between_sample_dict.get(sample_pair_key,-1):
                    between_sample_dict[sample_pair_key] = percent_similarity
                    #Check if the percent simliarites are greater than the threshold
                    if percent_similarity >= linkage_threshold:
                        #create haplotype pair alphabetically
                        haplotype_key = (sample_name,hit_sample
                            ) if sample_name < hit_sample else (hit_sample,sample_name)
                        #store linkage pair
                        haplotype_linkage_dict[sample_pair_key] = (haplotype_key,percent_similarity)
        
    #return tuple of each of the dictoinaries in a tuple for the multiprocessing worker
    return (within_sample_dict, between_sample_dict, haplotype_linkage_dict)

def compare_and_compute_similarities(hash_table,
                                     samples_dir_to_compare,
                                     kmer_size,
                                     linkage_threshold):
    """
    Compare samples in the input directory to the hash table and compute the within and between sample percent simliarties. Also look for potential linkages.
    
    Args:
        hash_table (FullHashTable): The hash table to search for k-mers.
        samples_dir_to_compare (str): The directory containing the fasta files to compare.
        kmer_size (int): The size of the k-mers to be used in the analysis.
        linkage_threshold (int): The percent similarity threshold to store specific haplotype linkages
    Returns:
        tuple: A tuple containing three dictionaries:
            - within_sample_dict (dict): A dictionary with keys as tuples of sample pairs and values as their percent similarity for within-sample comparisons.
            - between_sample_dict (dict): A dictionary with keys as tuples of sample pairs and values as their percent similarity for between-sample comparisons.
            - haplotype_linkage_dict (dict): A dictionary with keys as tuples of sample pairs and values as tuples containing the best haplotype pair and their percent similarity for potentially linked samples.
    """

    #get the total unique kmers for all haplotypes in the hash table
    id_counts = hash_table.get_sample_id_counts()

    #collect all fasta files from the input directory
    fasta_files = [ f for f in os.listdir(samples_dir_to_compare)
                   if f.endswith(".fasta") or f.endswith(".fa") ]

    #prepare an argument list containing a tuple per fasta file
    args_list = [( fasta,
                samples_dir_to_compare,
                kmer_size,
                linkage_threshold)
                for fasta in fasta_files]
    #initiailzie final dictionaries to store results
    within_sample_dict = {}
    between_sample_dict = {}
    haplotype_linkage_dict = {}
    # Run processing in parallel (one FASTA file per worker process)
    with Pool(cpu_count(),
            initializer=init_worker,        # initialize global hash table + ID counts once per worker
            initargs=(hash_table, id_counts)) as pool:

        # Iterate over completed FASTA jobs as they finish
        # imap_unordered yields results immediately when a worker completes
        # tqdm wraps the iterator to display a progress bar
        for w, b, h in tqdm(
                pool.imap_unordered(process_single_fasta, args_list),
                total=len(args_list),          # total number of FASTA files to process
                desc="Processing FASTAs",     # progress bar label
                unit="fasta"):                # unit displayed in progress bar

            ### MERGE WITHIN-SAMPLE RESULTS ###
            # w is the within_sample_dict returned from one worker
            # safe to update directly because keys are haplotype pairs
            within_sample_dict.update(w)

            ### MERGE BETWEEN-SAMPLE RESULTS ###
            # b is the between_sample_dict returned from one worker
            # we keep only the highest percent similarity per sample pair
            for k, v in b.items():
                # if this similarity is higher than what we’ve seen before, replace it
                if v > between_sample_dict.get(k, -1):
                    between_sample_dict[k] = v

            ### MERGE HAPLOTYPE LINKAGE RESULTS ###
            # h stores:
            # {sample_pair : ((hap1, hap2), percent_similarity)}
            # we again keep only the highest percent similarity per sample pair
            for k, v in h.items():
                # v[1] is the percent similarity
                # compare against stored similarity (default = -1 if not yet present)
                if v[1] >= haplotype_linkage_dict.get(k, (None, -1))[1]:
                    haplotype_linkage_dict[k] = v

    #conver from :
    #{sample_pair : ((haplotype_1,hapltype_2),percent_similarity)}
    #to: 
    #{(haplotype_1,hapltype_2):percent_similarity}
    haplotype_linkage_dict = {
        v[0]: v[1]
        for v in haplotype_linkage_dict.values()
    }
    return within_sample_dict, between_sample_dict, haplotype_linkage_dict


#function to run the hash table initialization and key insertion
def update_table(previous_hash_table, samples_dir_to_add, new_hash_table_name, save_table,kmer_size):
    """
    Updates an existing hash table with new samples and optionally saves the updated table.

    Args:
        previous_hash_table (str): Path to the file containing the previous hash table.
        samples_dir_to_add (str): Directory containing new samples to add to the hash table.
        new_hash_table_name (str): Name for the new hash table file to be saved.
        save_table (bool): Flag indicating whether to save the updated hash table to a file.

    Returns:
        FullHashTable: The updated hash table object.
    """
    # Load the previous hash table
    hash_table = FullHashTable.load_table(previous_hash_table)
    # Iterate over the samples in the directory
    hash_table = iterate_over_samples(samples_dir_to_add, hash_table,kmer_size)
    # Save the hash table to a file
    if save_table:
        hash_table.save_table(f"{new_hash_table_name}.pkl")
    return hash_table
def create_new_table(samples_dir, new_hash_table_name, save_table, kmer_size):
    """
    Create a new hash table, populate it with data from samples, and optionally save it to a file.

    Args:
        samples_dir (str): The directory containing sample data to be hashed.
        new_hash_table_name (str): The name to be used when saving the hash table.
        save_table (bool): A flag indicating whether to save the hash table to a file.
        kmer_size (int): The size of the k-mers to be used in the analysis.

    Returns:
        FullHashTable: The populated hash table.
    """
    #create a new hash table
    hash_table = FullHashTable()
    #iterate over the samples in the directory
    hash_table = iterate_over_samples(samples_dir, hash_table, kmer_size)
    #save the hashtable to a file
    if save_table:
        hash_table.save_table(f"{new_hash_table_name}.pkl")
    return hash_table
def main(argv=None): 
    #parse command line arguments
    parser = argparse.ArgumentParser(description='Hash table and kmer-ization implementation for HCV clustering')
    parser.add_argument('--mode', type=str, choices=['initialize', 'update', 'load'], required=True, help='Mode to perform: initialize, update, or load')
    parser.add_argument('--samples_dir_to_add', type=str, help='Directory of samples to add to the hash table')
    parser.add_argument('--samples_dir_to_compare', type=str, help='Directory of samples to compare to the hash table')
    parser.add_argument('--previous_hash_table', type=str, help="Path to a previous hash table to load in and use")
    parser.add_argument('--new_hash_table_name', type=str, default='hcv_hash_table', help='Name of the final outputted Hash table. If a hash table is being updated, the updated table will be saved to the new name')
    parser.add_argument('--kmer_size', type=int, default=10, help='Size of the k-mer to be used in the analysis. Default is 10')
    parser.add_argument('--save_table', action='store_true', help='Save the hash table after creation')
    parser.add_argument('--linkage_threshold',type=int,default=50,help="Percent similarity threshold to store the specific haplotype linakges of potentially linked samples. Default:50" )
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args(argv)
    
    logging.info('Starting code\n')
    hash_table = None
    #print(the kmer size being used)
    logging.info(f"Kmer size being used is {args.kmer_size}")
    if args.mode == 'initialize':
        if args.samples_dir_to_add is None:
            logging.error("If --mode is 'initialize', --samples_dir_to_add (path) must be provided.")
            exit(1)
        print("Initialize mode selected: a new hash table will be created from scratch\n")
        hash_table = create_new_table(args.samples_dir_to_add, args.new_hash_table_name, args.save_table, args.kmer_size)
    elif args.mode == 'update':
        if args.samples_dir_to_add is None:
            logging.error("If --mode is 'update', --samples_dir_to_add (path) must be provided.")
            exit(1)
        if args.previous_hash_table is None:
            logging.error("If --mode is 'update', --previous_hash_table (path) must be provided.")
            exit(1)
        print("Update mode selected: a previous hash table will be loaded and updated with new samples\n")
        hash_table = update_table(args.previous_hash_table, args.samples_dir_to_add, args.new_hash_table_name, args.save_table, args.kmer_size)
    elif args.mode == 'load':
        if args.previous_hash_table is None:
            logging.error("If --mode is 'load', --previous_hash_table (path) must be provided.")
            exit(1)
        print("Load mode selected: a previous hash table will be loaded\n")
        hash_table = FullHashTable.load_table(args.previous_hash_table)
    
    #if samples_dir_to_compare is provided, compare the samples 
    if args.samples_dir_to_compare is not None:
        print("Comparing samples in specified directory to hash table")
        logging.info(f"Using a linkage_threshold of : {args.linkage_threshold}")
        within_sample_dict, between_sample_dict, haplotype_linkage_dict = compare_and_compute_similarities(hash_table, args.samples_dir_to_compare, args.kmer_size, args.linkage_threshold)
        #generate the dataframes
        between_samples_df = pd.DataFrame(
            [(primary_sample, hit_sample, percent_similarity) for (primary_sample, hit_sample), percent_similarity in between_sample_dict.items()],
                    columns=["sample_1", "sample_2", "percent_similarity"]
        )
        within_samples_df = pd.DataFrame(
            [(primary_sample, hit_sample, percent_similarity) for (primary_sample, hit_sample), percent_similarity in within_sample_dict.items()],
                    columns=["sample_1", "sample_2", "percent_similarity"]
        )
        haplotype_linkage_df = pd.DataFrame(
            [(sample_1, sample_2, percent_similarity) for (sample_1,sample_2),percent_similarity in haplotype_linkage_dict.items()],
                    columns=["sample_1", "sample_2", "percent_similarity"]
        )
        #sort values by percent similarity for the three dataframes
        between_samples_df = between_samples_df.sort_values(by=["percent_similarity", "sample_1","sample_2"], ascending=False)
        within_samples_df = within_samples_df.sort_values(by=["percent_similarity", "sample_1","sample_2"], ascending=False)
        haplotype_linkage_df = haplotype_linkage_df.sort_values(by=["percent_similarity", "sample_1","sample_2"], ascending=False)
        #save the two dataframes
        between_samples_df.to_csv('between_sample_percent_similarities.csv',index=False)
        within_samples_df.to_csv('within_sample_percent_similarities.csv',index=False)
        haplotype_linkage_df.to_csv('haplotype_linkages.csv',index=False)
    logging.info('Code finished running')
#write main function to test the class
if __name__ == "__main__":
    main()
