#!/usr/bin/env python
# -*- coding: utf-8 -*-
# hcv_hasher.py
# A script to implement a hash table for storing k-mer information sample IDs
import mmh3  # MurmurHash3 for hashing
import pandas as pd 
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import os 
import random
import pickle
import argparse
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class FullHashTable:
    def __init__(self, size):
        self.size = size
        self.table = np.full(size, None, dtype=object)  # Store (key, tax_id_set) tuples

    def insert(self, key, tax_id):
        """Insert a key-value pair using full hash storage."""
        full_hash = mmh3.hash(key, signed=False)  # Compute full 32-bit MurmurHash3
        index = full_hash % self.size  # Compute index using modulo operation

        # Linear probing for collision resolution
        while self.table[index] is not None:
            stored_key, stored_tax_ids = self.table[index]
            if stored_key == key:
                stored_tax_ids.add(tax_id)  # Add tax_id if key already exists
                return
            index = (index + 1) % self.size  # Move to the next index

        # Store key and associated tax_id set
        self.table[index] = (key, {tax_id})  

    def search(self, key):
        """Search for a key and return the associated taxonomic IDs."""
        full_hash = mmh3.hash(key, signed=False)
        index = full_hash % self.size

        # Linear probing for search
        while self.table[index] is not None:
            stored_key, stored_tax_ids = self.table[index]
            if stored_key == key:  # Fix: Compare key directly
                return stored_tax_ids  # Return the set of taxonomic IDs
            index = (index + 1) % self.size
        return None  # Key not found

    def resize(self, new_size):
        """Resize the hash table and rehash all elements."""
        temp_hash_table = FullHashTable(new_size)  # Create a new hash table
        for entry in self.table:
            if entry is not None:
                key, tax_ids = entry
                for tax_id in tax_ids:  
                    temp_hash_table.insert(key, tax_id)  # Reinsert correctly
        # Update to new table
        self.size = temp_hash_table.size
        self.table = temp_hash_table.table

    def count_filled(self):
        """Count the number of filled and empty slots in the hash table."""
        filled = sum(1 for entry in self.table if entry is not None)
        empty = self.size - filled
        #print(f"Slots filled: {filled}, empty slots: {empty}")
        return filled, empty

    def save_table(self, filename):
        """Save the hash table to a file using pickle."""
        with open(filename, "wb") as f:
            pickle.dump((self.size, self.table), f)
        print(f"Hash table saved to {filename}")

    @staticmethod
    def load_table(filename):
        """Load a hash table from a file."""
        with open(filename, "rb") as f:
            size, table = pickle.load(f)
        new_ht = FullHashTable(size)
        new_ht.table = table
        print(f"Hash table loaded from {filename}")
        return new_ht

    def get_tax_id_counts(self):
        """Count occurrences of each taxonomic ID in the hash table."""
        tax_id_counts = {}
        for entry in self.table:
            if entry is not None:
                _, tax_ids = entry
                for tax_id in tax_ids:
                    tax_id_counts[tax_id] = tax_id_counts.get(tax_id, 0) + 1
        return tax_id_counts

    def get_tax_id_counts_df(self):
        """Return a DataFrame with tax_id and their respective counts."""
        tax_id_counts = self.get_tax_id_counts()
        return pd.DataFrame(list(tax_id_counts.items()), columns=["tax_id", "counts"])
    def display_table(self):
        print(self.table)
        
# Produce a list of k-mer sets for each sequence in the input list
def produce_kmers(seq_list, k_mer=26):
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

def iterate_over_samples(samples_dir, hash_table):
    for filename in os.listdir(samples_dir):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            filepath = os.path.join(samples_dir, filename)
            seq_list = [str(record.seq) for record in SeqIO.parse(filepath, "fasta")]
            print(f"Processing {filename} with {len(seq_list)} sequences.")
            #produce the kmer sets
            k_mer_sets, num_sequences = produce_kmers(seq_list)
            #set the unique haplotype number
            haplotype_number = 0
            #now iterate through each set of kmers
            for kmer_set in k_mer_sets:
                #create the sample name 
                sample_name = f"{filename}_{haplotype_number}"
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
def compare_samples(samples_dir_to_compare, hash_table):
    pass
#function to run the hash table initialization and key insertion
def update_table(previous_hash_table, samples_dir_to_add, hash_table_name, save_table):
    # Load the previous hash table
    hash_table = FullHashTable.load_table(previous_hash_table)
    print(f"Size of loaded table{hash_table.size}")
    # Iterate over the samples in the directory
    hash_table = iterate_over_samples(samples_dir_to_add, hash_table)
    # Save the hash table to a file
    if save_table:
        hash_table.save_table(f"{hash_table_name}.pkl")
    return hash_table
def create_new_table(table_size, samples_dir, table_name,save_table):
    #create a new hash table
    hash_table = FullHashTable(table_size)
    #iterate over the samples in the directory
    hash_table = iterate_over_samples(samples_dir, hash_table)
    #save the hashtable to a file
    if save_table:
        hash_table.save_table(f"{table_name}.pkl")
    return hash_table
def main(argv=None):
    logging.info('Starting code\n')
    #parse command line arguments
    parser = argparse.ArgumentParser(description='Hash table and kmer-ization implementation for HCV clustering')
    parser.add_argument('--table_size', type=int, default=500000, help='Size of the hash table to initialize with')
    parser.add_argument('--mode', type=str, choices=['new', 'update', 'load'], required=True, help='Mode to perform: new, update, or load')
    parser.add_argument('--samples_dir_to_add', type=str, help='Directory of samples to add to the hash table')
    parser.add_argument('--samples_dir_to_compare', type=str, help='Directory of samples to compare to the hash table')
    parser.add_argument('--previous_hash_table', type=str, help="Path to a previous hash table to load in and use")
    parser.add_argument('--hash_table_name', type=str, default='hcv_hash_table', help='Name of the final outputted Hash table. If a hash table is being updated, the updated table will be saved to the new name')
    parser.add_argument('--kmer_size', type=int, default=26, help='Size of the k-mer to be used in the analysis')
    parser.add_argument('--save_table', action='store_true', help='Save the hash table after creation')
    args = parser.parse_args(argv)
    
    hash_table = None
    
    if args.mode == 'new':
        if args.samples_dir_to_add is None:
            logging.error("If --mode is 'new', --samples_dir (path) must be provided.")
            exit(1)
        hash_table = create_new_table(args.table_size, args.samples_dir_to_add, args.hash_table_name, args.save_table)
    elif args.mode == 'update':
        if args.samples_dir_to_add is None:
            logging.error("If --mode is 'update', --samples_dir (path) must be provided.")
            exit(1)
        if args.previous_hash_table is None:
            logging.error("If --mode is 'update', --previous_hash_table (path) must be provided.")
            exit(1)
        hash_table = update_table(args.previous_hash_table, args.samples_dir_to_add, args.hash_table_name, args.save_table)
    elif args.mode == 'load':
        if args.previous_hash_table is None:
            logging.error("If --mode is 'load', --previous_hash_table (path) must be provided.")
            exit(1)
        hash_table = FullHashTable.load_table(args.previous_hash_table)
    
    # if args.samples_dir is not None and args.mode != 'load':
    #     compare_samples(args.samples_dir, hash_table)
    
    logging.info('Code finished running')
#write main function to test the class
if __name__ == "__main__":
    main()