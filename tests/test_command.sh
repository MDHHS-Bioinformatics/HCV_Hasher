#!/bin/bash

#a simple test command to get used to running HCV Hasher
#compare outputs of this command to the results in : tests/data/results
#to ensure reproducibility
python hcv_hasher.py --mode initialize \
    --samples_dir_to_add tests/data/fastas/ \
    --samples_dir_to_compare tests/data/fastas/