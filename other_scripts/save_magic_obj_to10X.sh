#!/bin/bash -l
#$ -N svae_parse_magic_to_10X
# specify queue
#$ -q all.q
# start in the current working directory
#$ -cwd
# email notification
#$ -M rwo8@cdc.gov
# email if job aborts
#$ -m a
##$ -t 1-3
# N 12
# mem=50G

conda activate python3

Rscript save_magic_parse_into_10X.r

conda deactivate