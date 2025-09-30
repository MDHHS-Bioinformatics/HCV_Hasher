#!/usr/bin/env python
# -*- coding: utf-8 -*-
# hcv_hasher.py
# A script to implement a hash table for storing k-mer information sample IDs

__version__ = "1.0.2"

import mmh3  # MurmurHash3 for hashing
import pandas as pd 
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import os 
import pickle
import argparse
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class FullHashTable:
    """
    A class to represent a hash table for storing k-mer information with associated sample IDs.
    """

    def __init__(self, size):
        """
        Initialize the hash table with a given size.

        Parameters:
        size (int): The size of the hash table.
        """
        self.size = size
        self.table = np.full(size, None, dtype=object)  # Store (key, sample_id_set) tuples

    def insert(self, key, sample_id):
        """
        Insert a key-value pair into the hash table.

        Parameters:
        key (str): The k-mer to be inserted.
        sample_id (str): The sample ID associated with the k-mer.
        """
        full_hash = mmh3.hash(key, signed=False)  # Compute full 32-bit MurmurHash3
        index = full_hash % self.size  # Compute index using modulo operation

        # Linear probing for collision resolution
        while self.table[index] is not None:
            stored_key, stored_sample_ids = self.table[index]
            if stored_key == key:
                stored_sample_ids.add(sample_id)  # Add sample_id if key already exists
                return
            index = (index + 1) % self.size  # Move to the next index

        # Store key and associated sample_id set
        self.table[index] = (key, {sample_id})  

    def search(self, key):
        """
        Search for a key in the hash table and return the associated sample IDs.

        Parameters:
        key (str): The k-mer to search for.

        Returns:
        set: A set of sample IDs associated with the k-mer, or None if the key is not found.
        """
        full_hash = mmh3.hash(key, signed=False)
        index = full_hash % self.size

        # Linear probing for search
        while self.table[index] is not None:
            stored_key, stored_sample_ids = self.table[index]
            if stored_key == key:  # Fix: Compare key directly
                return stored_sample_ids  # Return the set of taxonomic IDs
            index = (index + 1) % self.size
        return None  # Key not found

    def resize(self, new_size):
        """
        Resize the hash table and rehash all elements.

        Parameters:
        new_size (int): The new size of the hash table.
        """
        temp_hash_table = FullHashTable(new_size)  # Create a new hash table
        for entry in self.table:
            if entry is not None:
                key, sample_ids = entry
                for sample_id in sample_ids:  
                    temp_hash_table.insert(key, sample_id)  # Reinsert correctly
        # Update to new table
        self.size = temp_hash_table.size
        self.table = temp_hash_table.table

    def count_filled(self):
        """
        Count the number of filled and empty slots in the hash table.

        Returns:
        tuple: A tuple containing the number of filled slots and empty slots.
        """
        filled = sum(1 for entry in self.table if entry is not None)
        empty = self.size - filled
        return filled, empty

    def save_table(self, filename):
        """
        Save the hash table to a file using pickle.

        Parameters:
        filename (str): The name of the file to save the hash table to.
        """
        with open(filename, "wb") as f:
            pickle.dump((self.size, self.table), f)
        print(f"Hash table saved to {filename}")

    @staticmethod
    def load_table(filename):
        """
        Load a hash table from a file.

        Parameters:
        filename (str): The name of the file to load the hash table from.

        Returns:
        FullHashTable: The loaded hash table.
        """
        with open(filename, "rb") as f:
            size, table = pickle.load(f)
        new_ht = FullHashTable(size)
        new_ht.table = table
        print(f"Hash table loaded from {filename}")
        return new_ht

    def get_sample_id_counts(self):
        """
        Count occurrences of each sample ID in the hash table.

        Returns:
        dict: A dictionary with sample IDs as keys and their counts as values.
        """
        sample_id_counts = {}
        for entry in self.table:
            if entry is not None:
                _, sample_ids = entry
                for sample_id in sample_ids:
                    sample_id_counts[sample_id] = sample_id_counts.get(sample_id, 0) + 1
        return sample_id_counts

    def get_sample_id_counts_df(self):
        """
        Return a DataFrame with sample IDs and their respective counts.

        Returns:
        pandas.DataFrame: A DataFrame with columns "sample_id" and "counts".
        """
        sample_id_counts = self.get_sample_id_counts()
        return pd.DataFrame(list(sample_id_counts.items()), columns=["sample_id", "counts"])

    def display_table(self):
        """
        Display the contents of the hash table.
        """
        print(self.table)
        
# Produce a list of k-mer sets for each sequence in the input list
def produce_kmers(seq_list, k_mer=25):
    """
    Generate k-mers for each sequence in the provided list.

    Args:
        seq_list (list of str): List of sequences to process.
        k_mer (int, optional): Length of the k-mers to generate. Default is 25.

    Returns:
        tuple: A tuple containing:
            - list of set: A list where each element is a set of k-mers for the corresponding sequence.
            - int: The total number of sequences processed.
    """
    k_mer_sets = []  # List to hold sets of k-mers for each sequence
    # Iterate over each sequence in the list
    for seq in seq_list:
        k_mer_set = set()  # Create a new set for the current sequence
        for i in range(len(seq) - k_mer + 1):
            k_mer_set.add(seq[i:i + k_mer])
        k_mer_sets.append(k_mer_set)  # Append the set to the list
    #print(f"Total sequences processed: {len(k_mer_sets)}")
    #return the total sequences as well
    return k_mer_sets, len(k_mer_sets)

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
            fasta_name = filename.split('.')[0]
            seq_list = [str(record.seq) for record in SeqIO.parse(filepath, "fasta")]
            print(f"Processing {fasta_name} with {len(seq_list)} sequences.")
            #produce the kmer sets
            k_mer_sets, num_sequences = produce_kmers(seq_list,k_mer=kmer_size)
            #set the unique haplotype number
            haplotype_number = 0
            #now iterate through each set of kmers
            for kmer_set in k_mer_sets:
                #create the sample name 
                sample_name = f"{fasta_name}_{haplotype_number}"
                #add each kmer to the hash table
                for kmer in kmer_set:
                    hash_table.insert(kmer, sample_name)
                #increment the haplotype number
                haplotype_number+=1
            #after every sample added check the load size
            filled,empty_slots = hash_table.count_filled()
            print(f"filled: {filled}, empty_slots: {empty_slots}")
            load_factor = filled/(filled + empty_slots)
            print(f"Load factor: {load_factor}")
            #if the filled slots are greater than 65% of the total size, resize the table
            if load_factor > 0.65:
                print("Resizing hash table...")
                hash_table.resize(hash_table.size * 2)
                print(f"Table resized to {hash_table.size}")
    print("\nAll samples processed and added to hash table")
    return hash_table
#Function to store only the sample sample comarisions in a condensed dictionary
def get_within_sample_values(threshold_dict):
    """
    Extracts and returns a dictionary of threshold values for sample pairs that belong to the same base sample.

    This function iterates over a dictionary where keys are tuples of sample names and values are percent similarity values.
    For each pair, it checks if the base names of the samples (defined as the part of the sample name before the last underscore)
    are the same but the full sample names are different. If so, it adds the sorted pair and its value to the result dictionary,
    ensuring that each pair is only included once (regardless of order).

    Args:
        threshold_dict (dict): A dictionary with keys as tuples of sample names (str, str) and values as percent similarities.

    Returns:
        dict: A dictionary containing only those pairs where the base sample names match but the full sample names differ,
              with keys as sorted tuples of sample names and values as the corresponding threshold values.
    """
    #initialize a new dictionary to store our values
    within_sample_threshold_dict = {}
    #iterate over the the threshold dictionary
    for (sample_1, sample_2), value in threshold_dict.items():
        # Get the base names of the samples (everything before the last underscore)
        sample_1_base = sample_1.rsplit('_', 1)[0]
        sample_2_base = sample_2.rsplit('_', 1)[0]
        #Check if the base names are the same
        if sample_1_base == sample_2_base:
            #if the sample names are the same, skip this iteration
            if sample_1 == sample_2:
                continue
            #sort keys in more efficient manner
            if sample_1 < sample_2:
                new_key = (sample_1, sample_2)
            else:
                new_key = (sample_2, sample_1)
            #store the hapltoype-hapltype comparison for each sample
            within_sample_threshold_dict[new_key] = value
    return within_sample_threshold_dict
    
    
#Function to store only unique sample comparitons and keep the most simliar pairs between samples
def get_between_sample_values(threshold_dict,threshold=50):
    """
    Processes a dictionary of sample pair similarity values and returns a new dictionary
    containing only the highest similarity value for each unique pair of base sample names,
    excluding within-sample comparisons.

    Args:
        threshold_dict (dict): A dictionary where keys are tuples of sample names (str, str)
            and values are similarity scores (numeric). 

    Returns:
        dict: A dictionary with keys as tuples of base sample names (str, str), sorted
            alphabetically, and values as the highest similarity score observed between
            any pair of samples from those base names.
        dict: A dictionary with keys as tuples of original sample names (str, str) that meet or exceed
            the specified threshold, and values as their corresponding similarity scores.

    Notes:
        - Only between-sample comparisons are included (i.e., pairs with different base names).
        - For each unique pair of base sample names, only the maximum similarity value is kept.
        - The base sample name is defined as the portion of the sample name before the last underscore.
    """
    #initialize a new dictionary to store our values 
    between_sample_dict = {}
    haplotype_linkage_dict = {}
    for (sample_1,sample_2),value in threshold_dict.items():
        #clean the sample names to just use the base name (everything before the last underscore
        sample_1_base = sample_1.rsplit('_', 1)[0]
        sample_2_base = sample_2.rsplit('_', 1)[0]
        #only do comparisons between samples not within samples
        if sample_1_base == sample_2_base:
            continue # skipping same sample comparison 
        
        #sort keys in more efficient manner
        if sample_1_base < sample_2_base:
            new_key = (sample_1_base, sample_2_base)
            inner_key = (sample_1, sample_2)
        else:
            new_key = (sample_2_base, sample_1_base)
            inner_key = (sample_2, sample_1)
            
        #keep max similarity value for each sample pair
        if value > between_sample_dict.get(new_key, -float('inf')):
            between_sample_dict[new_key] = value
            if value >= threshold:
                inner_key+= (value,)
                haplotype_linkage_dict[new_key] = inner_key
                
    return between_sample_dict, haplotype_linkage_dict  

def get_similarity_values(sample_search_results, sample_kmer_counts, hash_table):
    """
    Calculate the similarity values between samples based on k-mer counts and a given threshold.

    Args:
        sample_search_results (dict): A dictionary where keys are primary sample IDs and values are dictionaries 
                                          of hit sample IDs and their respective k-mer hit counts.
        sample_kmer_counts (dict): A dictionary where keys are sample IDs and values are their respective total k-mer counts.
        hash_table (object): An object that contains the method get_sample_id_counts() which returns a dictionary of sample IDs 
                                 and their respective total k-mer counts.
        threshold (int, optional): The minimum percentage similarity required to include the sample pair in the result. 
                                       Defaults to 50.

    Returns:
        dict: A dictionary where keys are tuples of (primary_sample, hit_sample) and values are the percentage similarity 
                  between the primary sample and the hit sample, only including pairs that meet or exceed the threshold.
    """
    id_counts = hash_table.get_sample_id_counts()
    #create a new threshold dictionary to store samples who meet the threshold 
    threshold_dict = {}  
    for primary_sample, search_results in sample_search_results.items():
        #get the kmer total for the primary sample
        primary_sample_kmer_count = sample_kmer_counts[primary_sample]
        #iterate over all the hits that we got in the hash table 
        for hit_sample, hit_counts in search_results.items():
            #get the total counts for sample hit from the hash table
            hit_sample_kmer_count = id_counts[hit_sample]
            #calculate the total unique kmer counts between our primary and hit sample 
            total_unique_counts = (primary_sample_kmer_count + hit_sample_kmer_count) - hit_counts
            #calculate the percent simliitariy and round to three decimal places
            percent_similarity = round(((hit_counts / total_unique_counts) * 100),3)
            #if the percent similarity is greater than the threshold, store it in a new dictionary
            #if percent_similarity >= threshold:
            threshold_dict[(primary_sample, hit_sample)] = percent_similarity
    return threshold_dict
def compare_samples(hash_table,samples_dir_to_compare,kmer_size):
    """
    Compare samples in a given directory against a hash table of kmers.
    Args:
        hash_table (object): An object that supports a `search` method for kmers.
        samples_dir_to_compare (str): Path to the directory containing sample fasta files.
    Returns:
        tuple: A tuple containing two dictionaries:
            - sample_search_results (dict): A dictionary where keys are sample names and values are dictionaries 
                of other samples and their respective kmer counts.
            - sample_kmer_counts (dict): A dictionary where keys are sample names and values are the counts of kmers 
                for each haplotype in the sample.
    Notes:
        - Only fasta files (with extensions .fasta or .fa) in the samples directory are processed.
        - The sample name is constructed from the fasta file name and haplotype number.
        - The function does not store search results for the current sample if it is found in the hash table.
        """
    #create dictionary to store what toher samples each kmer belonged too for each haplotype
    sample_search_results = {}
    #create a dictionary to store the kmer coutns for each haplotype for each sample
    sample_kmer_counts = {}
    #iterate over the samples directory 
    for fasta in os.listdir(samples_dir_to_compare):
        #only process fasta files
        if fasta.endswith(".fasta") or fasta.endswith(".fa"):
            fasta_name = fasta.split('.')[0]
            #get the genotype 
            genotype = fasta_name.split('-')[-1]
            fasta_path = os.path.join(samples_dir_to_compare, fasta)
            seq_list = [str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")]
            #get the kmers for the sample
            k_mer_sets, num_sequences = produce_kmers(seq_list,k_mer=kmer_size)
            #set the unique haplotype number
            haplotype_number = 0
            #now iterate through each set of kmers
            for kmer_set in k_mer_sets:
                #create a temp dictionary to store the search results
                temp_dict = {}
                kmer_count=0
                #create the sample name 
                sample_name = f"{fasta_name}_{haplotype_number}"
                #search for each kmer in the hash table
                for kmer in kmer_set:
                    #increase the kmer count
                    kmer_count+=1
                    result_set = hash_table.search(kmer)
                    #print(result_set)
                    #iterate through each sampe name in the results set 
                    if result_set is not None:
                        for result in result_set:
                            #if fasta_name not in result: #this will prevent the current sample from being stored 
                            if result in temp_dict:
                                temp_dict[result] += 1
                            else:
                                temp_dict[result] = 1
                #update the respective dictionaries
                sample_kmer_counts[sample_name] = kmer_count
                #store the search results in the dictionary
                sample_search_results[sample_name] = temp_dict
                #increment the haplotype number     
                haplotype_number+=1
                    
    #return the dictionaries
    return sample_search_results, sample_kmer_counts
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
    print(f"Size of loaded table: {hash_table.size}")
    # Iterate over the samples in the directory
    hash_table = iterate_over_samples(samples_dir_to_add, hash_table,kmer_size)
    # Save the hash table to a file
    if save_table:
        hash_table.save_table(f"{new_hash_table_name}.pkl")
    return hash_table
def create_new_table(table_size, samples_dir, new_hash_table_name, save_table, kmer_size):
    """
    Create a new hash table, populate it with data from samples, and optionally save it to a file.

    Args:
        table_size (int): The size of the hash table to be created.
        samples_dir (str): The directory containing sample data to be hashed.
        new_hash_table_name (str): The name to be used when saving the hash table.
        save_table (bool): A flag indicating whether to save the hash table to a file.

    Returns:
        FullHashTable: The populated hash table.
    """
    #create a new hash table
    hash_table = FullHashTable(table_size)
    #iterate over the samples in the directory
    hash_table = iterate_over_samples(samples_dir, hash_table, kmer_size)
    #save the hashtable to a file
    if save_table:
        hash_table.save_table(f"{new_hash_table_name}.pkl")
    return hash_table
def main(argv=None): 
    #parse command line arguments
    parser = argparse.ArgumentParser(description='Hash table and kmer-ization implementation for HCV clustering')
    parser.add_argument('--table_size', type=int, default=500000, help='Size of the hash table to initialize with')
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
        hash_table = create_new_table(args.table_size, args.samples_dir_to_add, args.new_hash_table_name, args.save_table, args.kmer_size)
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
        sample_search_results, sample_kmer_counts = compare_samples(hash_table, args.samples_dir_to_compare, args.kmer_size)
        #get the similarity values
        similarity_dict = get_similarity_values(sample_search_results,sample_kmer_counts,hash_table)
        #get the within sample similarity values
        within_sample_dict = get_within_sample_values(similarity_dict)
        #get the between sample similarity values 
        between_sample_dict,haplotype_linkage_dict = get_between_sample_values(similarity_dict,args.linkage_threshold)
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
            [(sample_1_base, sample_2_base, sample_1, sample_2, percent_similarity) for (sample_1_base, sample_2_base), (sample_1,sample_2,percent_similarity) in haplotype_linkage_dict.items()],
                    columns=["sample_1_base", "sample_2_base", "sample_1", "sample_2", "percent_similarity"]
        )
        #sort values by percent similarity for the three dataframes
        between_samples_df = between_samples_df.sort_values(by=["percent_similarity", "sample_1","sample_2"], ascending=False)
        within_samples_df = within_samples_df.sort_values(by=["percent_similarity", "sample_1","sample_2"], ascending=False)
        haplotype_linkage_df = haplotype_linkage_df.sort_values(by=["percent_similarity", "sample_1","sample_2"], ascending=False)
        #drop sample_1_base and sample_2_base columns
        haplotype_linkage_df = haplotype_linkage_df.drop(columns=["sample_1_base", "sample_2_base"])
        #save the two dataframes
        between_samples_df.to_csv('between_sample_percent_similarities.csv',index=False)
        within_samples_df.to_csv('within_sample_percent_similarities.csv',index=False)
        haplotype_linkage_df.to_csv('haplotype_linkages.csv',index=False)
    logging.info('Code finished running')
#write main function to test the class
if __name__ == "__main__":
    main()