library(Rmagic, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(readr, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(dplyr, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(tibble)
library(Seurat, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(Matrix)
library(ggplot2, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(glmGamPoi, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(future, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(SeuratDisk, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(HGNChelper, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(DropletUtils, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")


fig_path <- file.path("/scicomp/home-pure/rwo8/workdir/single_cell/parse/output/figures")

table_path <- file.path("/scicomp/home-pure/rwo8/workdir/single_cell/parse/output/tables")

marker_path <- "/scicomp/home-pure/rwo8/workdir/single_cell/cell_markers"

RDS_path <- "/scicomp/home-pure/rwo8/workdir/single_cell/RDS"

thread=8

###########################################################################################################
#
# 1) read from parse dataset
#
###########################################################################################################
print("######################STEP1: Read data into Seurat Object#########################################")

rds_files <- list.files(RDS_path, full.names = TRUE, pattern="MON*")

sample_objs <- lapply(rds_files, function(x){
  y <- LoadH5Seurat(x, assay="RNA", reduction = FALSE, graphs=FALSE)
  ycol <- ncol(y@meta.data)
  y@meta.data$is.doublet <- y@meta.data[,ycol]
  y@meta.data <- y@meta.data[,-c(ycol,ycol-1, ycol-2, ycol-3, ycol-4,ycol-5)]
  y@meta.data$seq_sample <- basename(x)
  rna_quantile <- quantile(y@meta.data$nCount_RNA, seq(0,1,0.01))
  mt_quantile <- quantile(y@meta.data$percent.mt, seq(0,1,0.01))
  print(basename(x))
  print(paste0("number of RNA count minimum: ", as.numeric(rna_quantile[2][[1]]), "; ", "number of RNA count maximum: ", as.numeric(rna_quantile[length(rna_quantile)-1][[1]])))
  print(paste0("number of maximum mt percentage: ", as.numeric(mt_quantile[length(mt_quantile)-1][[1]])))
  parse_obj <- subset(y, subset = percent.mt < as.numeric(mt_quantile[length(mt_quantile)-1][[1]]) & nCount_RNA > as.numeric(rna_quantile[2][[1]]) & nCount_RNA < as.numeric(rna_quantile[length(rna_quantile)-1][[1]]))
  y
})

combined_obj <- Reduce(function(x,y) merge(x,y, merge.data=TRUE), sample_objs)

singlet_cells <- rownames(combined_obj@meta.data[combined_obj@meta.data[,"is.doublet"] == "Singlet",])
parse_obj <- subset(combined_obj, cells =singlet_cells)

metadata <- parse_obj@meta.data
metadata.1 <- metadata %>% mutate(cellID = mapply(function(x,y){
  ifelse(x == "MON14002A1_6000Cells.h5seurat", paste(y, "1", sep="_"),paste(y, "2", sep="_") )
},seq_sample,  bc_wells)) %>% mutate(Infection = sapply(sample, function(x){
  unlist(strsplit(x, "_", fixed=TRUE))[1]
}))
write_csv(metadata.1, "/scicomp/groups/OID/NCIRD/DVD/VGB/Norovirus/Rachel/single_cell/parseBioscience_pipeline/MON14002_01212022/partek_input_magic/cell_metadata.csv")
###################################################################################
#
# run MAGIC on raw dataset
#
###################################################################################
print("######################STEP4: Run magic #########################################")

# tv <- ifelse(infect == "NI",9 , ifelse(infect == "GII3", 10, ifelse(infect == "GII4", 11, 'auto')))
tv <- 'auto'
parse_obj <- magic(parse_obj, gene="all_genes", t = tv,  verbose=1, n.jobs=-1, seed = 123) # ALL combined cells t=11

write10xCounts(x = Matrix(parse_obj@assays$MAGIC_RNA@data, sparse = TRUE), path = "/scicomp/groups/OID/NCIRD/DVD/VGB/Norovirus/Rachel/single_cell/parseBioscience_pipeline/MON14002_01212022/partek_input_magic/")
