#!/bin/bash

#Run this before running anything else in the input directory
(ls | grep "R1_001.fastq") > input_R1_files.txt
(cut -d "_" -f 1 input_R1_files.txt) > sample_names.lst