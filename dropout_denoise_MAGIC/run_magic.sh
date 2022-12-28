#!/bin/bash -l
#$ -N report_denoise_raw_nig24
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

WORK="/scicomp/home-pure/rwo8/workdir/single_cell/script"

RDS_path="/scicomp/home-pure/rwo8/workdir/single_cell/parse/work_RDS"
Rscript $WORK/run_magic.r $RDS_path/raw_nig24_111922_normalized.RDS $RDS_path/raw_nig24_111922_magic.RDS


# Rscript $WORK/run_seurat_magic_combined.r 8 1.7 NI 20
# Rscript $WORK/run_seurat_noGastric.r 8 1.7 NI 20 4

# Rscript $WORK/run_seurat_magic_combined.r 8 1.7 GII3 30
# Rscript $WORK/run_seurat_noGastric.r 8 1.7 GII3 30 4

# Rscript $WORK/run_seurat_magic_combined.r 8 1.7 GII4 30
# Rscript $WORK/run_seurat_noGastric.r 8 1.7 GII4 30 4

################## plot ########################################################################

# Rscript $WORK/run.plot.umap.r

# Rscript $WORK/run.plot.tsne.r





