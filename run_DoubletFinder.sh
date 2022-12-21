#!/bin/bash -l
#$ -N run_doubletFinder_emory
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

################## run doublet finder ########################################################################
WORK="/scicomp/home-pure/rwo8/workdir/single_cell/script"

# ######## parse samples 
# PARSE="/scicomp/home-pure/rwo8/workdir/single_cell/MON14002/analysis"
# for dir in $PARSE/*00Cells;
# do
# sample=$(basename $dir)
# echo $sample
# Rscript $WORK/run_doublet_finder.r parse $dir/all-well/DGE_filtered/ $WORK/../RDS_2/$sample
# done

######### emory samples 
EMORY="/scicomp/home-pure/rwo8/workdir/single_cell/emory/Processed"
for dir in $EMORY/*;
do
sample=$(basename $dir)
echo $sample
Rscript $WORK/run_doublet_finder.r emory $dir/ $WORK/../RDS_2/$sample
done
