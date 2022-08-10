#!/bin/bash
#$ -M netid@nd.edu
#$ -m abe
#$ -r n
#$ -N BIP_2022_driver_jobOutput

# load the bio module for the servers
module load bio

# run the shell script for data prep and alignment
bash bip_2022_script.sh

# run the R script for quantifying transcript data
Rscript bip_2022_Rscript.R
