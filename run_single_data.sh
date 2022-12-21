#!/bin/bash -l
#$ -N emory_dataset
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

EMORY="/scicomp/home-pure/rwo8/workdir/single_cell/emory/Processed"
WORK="/scicomp/home-pure/rwo8/workdir/single_cell/script"

# for dir in $EMORY/*/;
# do
# sample=$(basename $dir)
# echo $sample
# Rscript $WORK/run_seurat_magic_emory.r $sample
# done

samples="NI GII3 GII4"
for sample in $samples;
do
echo $sample
Rscript $WORK/run_seurat_magic_parse.r $sample
done