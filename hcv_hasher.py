#!/usr/bin/env python
# -*- coding: utf-8 -*-
# hcv_hasher.py
# A script to implement a hash table for storing k-mer information sample IDs

__version__ = "1.0.3"

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
from collections import Counter 

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

def iter_kmers(seq, k):
    for i in range(len(seq) - k + 1):
        yield seq[i:i+k]

def compare_and_compute_similarities(hash_table,samples_dir_to_compare,kmer_size,linkage_threshold):
    """
    Memory-efficient combined comparison + similarity computation.
    """

    # Precompute total kmer counts from hash table
    id_counts = hash_table.get_sample_id_counts()

    #initiailize results dicitonaries
    within_sample_dict = {}
    between_sample_dict = {}
    haplotype_linkage_dict = {}

    #iterate over the fastas directory
    for fasta in os.listdir(samples_dir_to_compare):
        if not (fasta.endswith(".fasta") or fasta.endswith(".fa")):
            continue
        
        #construct the correct fasta path
        fasta_name = fasta.split('.')[0]
        fasta_path = os.path.join(samples_dir_to_compare, fasta)
        #set an intial haplotype number per sequence in the fasta file
        haplotype_number = 0
        #iterate through each sequence in the fasta file
        for record in SeqIO.parse(fasta_path, "fasta"):

            sample_name = f"{fasta_name}_{haplotype_number}"
            primary_base = sample_name.rsplit('_', 1)[0]

            #contains the number of shared kmers per individual sequences iterate through in the fasta
            temp_counter = Counter()
            kmer_count = 0
            seen_kmers = set() # Keep track of unique kmers to avoid double counting

            # Stream kmers (no sets stored)
            for kmer in iter_kmers(str(record.seq), kmer_size):
                if kmer in seen_kmers:
                    continue  # Skip already seen kmer
                seen_kmers.add(kmer)
                kmer_count += 1
                result_set = hash_table.search(kmer)
                #count the shared kmers
                if result_set:
                    for result in result_set:
                        temp_counter[result] += 1
            
            #compute percent similarities by iterating through each sample, and the shared number of kmers
            for hit_sample, hit_counts in temp_counter.items():
                #get the base sample name 
                hit_base = hit_sample.rsplit('_', 1)[0]
                hit_sample_kmer_count = id_counts.get(hit_sample,0)
                #calculate total unique counts 
                total_unique_counts = (
                    kmer_count + hit_sample_kmer_count - hit_counts
                )
                #calculate the percent similarity 
                percent_similarity = round((hit_counts / total_unique_counts) * 100, 3)
                
                ### WITHIN SAMPLE ###
                if primary_base == hit_base:
                    if sample_name == hit_sample:
                        continue #Skip self comparisons, they will obviously be 100% similiar
                    key = tuple(sorted([sample_name, hit_sample]))
                    within_sample_dict[key] = percent_similarity
                    
                ### BETWEEN SAMPLE ###
                else:
                    #alphabetical order the sample pair key
                    sample_pair_key =  (primary_base,hit_base 
                                        ) if primary_base < hit_base else (hit_base,primary_base)
                    
                    #check if this is the best similiarty so far for this pair 
                    if percent_similarity > between_sample_dict.get(sample_pair_key,-1):
                        between_sample_dict[sample_pair_key] = percent_similarity
                        #check if the percent similarity meets the linkage threshold
                        if percent_similarity >= linkage_threshold:
                            #store the best haplotype pair as the key, alphabetically ordered
                            haplotype_key = (sample_name, hit_sample
                                             ) if sample_name < hit_sample else (hit_sample, sample_name)
                            haplotype_linkage_dict[sample_pair_key] = ( haplotype_key, percent_similarity)
                            
                  
            #increase the haplotype number for the next sequence
            haplotype_number += 1
            #remove temporary variables to reset for next sequence
            del temp_counter
            del seen_kmers
            
    haplotype_linkage_dict = {v[0]: v[1] for v in haplotype_linkage_dict.values()}
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