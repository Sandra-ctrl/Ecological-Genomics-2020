#!/bin/bash

# we'll use this as a wrapper to run our different mapping scripts

# path to my repo:
myrepo="/users/s/n/snnadi/Ecological-Genomics-2020"

# My population:

mypop="KOS"

#Directory to our cleaned and paired reads:

input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${mypop}"

#Directory to store the output of our mapping
output="/data/project_data/RS_ExomeSeq/mapping"

#run mapping.sh

source ./mapping.sh

#run the post processing steps

source ./process_bam.sh 

