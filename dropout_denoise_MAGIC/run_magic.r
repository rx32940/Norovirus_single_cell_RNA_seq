library(Rmagic, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")

args = commandArgs(trailingOnly=TRUE)

input=args[1]
output=args[2]


# input="/scicomp/home-pure/rwo8/workdir/single_cell/parse/work_RDS/raw_nig24_111922_normalized.RDS"
# output="/scicomp/home-pure/rwo8/workdir/single_cell/parse/work_RDS/parse_harmony_magic_101722.RDS"

before_magic <- readRDS(input)

after_magic <- magic(before_magic, gene="all_genes", t = 'auto',  verbose=1, n.jobs=-1,seed = 47)

saveRDS(after_magic, file = output)

