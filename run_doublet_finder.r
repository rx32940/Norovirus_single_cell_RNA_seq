library(Rmagic, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(readr, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(dplyr, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(tibble)
library(Seurat, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(ggplot2, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(glmGamPoi, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(future, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")

args = commandArgs(trailingOnly=TRUE)

dataset <- args[1]
input <- args[2] # input path
output <- args[3]


# dataset <- "parse"
# input <- "/scicomp/home-pure/rwo8/workdir/single_cell/MON14002/analysis/MON14002A1_6000Cells/all-well/DGE_filtered/"
# output <- "/scicomp/home-pure/rwo8/workdir/single_cell/RDS/MON14002A1_6000Cells"
###########################################################################################################
#
# 1) read data
#
###########################################################################################################

read_parse <- function(input){
    DGE_folder <-input

    # gene expression matrix
    mat <- readMM(paste0(DGE_folder, "DGE.mtx"))

    # metadata and statistics of reads and genes of each well (cell)
    cell_meta <- read.delim(paste0(DGE_folder, "cell_metadata.csv"),
                    stringsAsFactor = FALSE, sep = ",")

    # reference gene annotation
    genes <- read.delim(paste0(DGE_folder, "all_genes.csv"),
                        stringsAsFactor = FALSE, sep = ",")
    
    # if gene does not have a gene name, replace with ensemble id
    genes[genes["gene_name"] == "",]$gene_name <- genes[genes["gene_name"] == "", "gene_id"]
    # if there is a duplication in cell bc id, make the id unique by adding _dup behind the id
    cell_meta$bc_wells <- make.unique(cell_meta$bc_wells, sep = "_dup")

    rownames(cell_meta) <- cell_meta$bc_wells # make the name of the cells (id), rownames
    genes$gene_name <- make.unique(genes$gene_name, sep = "_dup") # if duplicated genes, add _dup behind gene_name

    # Setting column and rownames to expression matrix
    colnames(mat) <- genes$gene_name # column of the matrix are gene names
    rownames(mat) <- rownames(cell_meta) # row of the matrix are cell (well) ids
    mat_t <- t(mat) # transpose the matrix so now rows are genes and columns are wells

    # Remove empty rownames, if they exist
    mat_t <- mat_t[(rowSums(mat_t) != 0),] # genes not expressed at any cell

    # if virus genes present in any of the cells
    rowSums(mat_t[rownames(mat_t) %in% c( "UHB38528", "UHB38529", "UHB38530","UHD96020","UHD96021", "UHD96022"),]) # quick check
    vir_genes <- mat_t[rownames(mat_t) %in% c( "UHB38528", "UHB38529", "UHB38530","UHD96020","UHD96021", "UHD96022"),]
    vir_df <- as.data.frame(t(vir_genes[,colSums(vir_genes) != 0])) %>% rownames_to_column()
    vir_df_cells <- cell_meta[rownames(cell_meta) %in% vir_df$rowname,] %>% rownames_to_column() %>% select(rowname,sample)
    vir_df_cells_sample <- left_join(vir_df,vir_df_cells, by="rowname" )
    write.csv(vir_df_cells_sample,file.path("/scicomp/home-pure/rwo8/workdir/single_cell/parse/output/tables", paste0(basename(output),"_virus_gene_expression_cells.csv" )))

    parse_obj <- CreateSeuratObject(mat_t, meta.data = cell_meta,  min.genes = 1, min.cells = 1)  

    parse_obj@meta.data$dataset <- factor(rep("Parse", nrow(parse_obj@meta.data)))

    parse_obj@meta.data$orig.ident <-factor(sapply(parse_obj@meta.data$sample, function(x){
    paste(unlist(strsplit(x, "_", fixed = TRUE))[1], "Parse", sep="_")
    }))
    Idents(parse_obj) <- parse_obj@meta.data$orig.ident

    parse_obj@meta.data$Sample <- factor(sapply(parse_obj@meta.data$sample, function(x){
    unlist(strsplit(x, "_", fixed = TRUE))[1]
    }))

    # percentage of mitohondria reads per cell
    parse_obj[["percent.mt"]] <- PercentageFeatureSet(parse_obj , pattern = "^MT-")
    # percentage of ribosomomal reads per cell
    parse_obj[["percent.rb"]] <- PercentageFeatureSet(parse_obj, pattern = "^RP[SL]")

    rna_quantile <- quantile(parse_obj@meta.data$nCount_RNA, seq(0,1,0.01))
    mt_quantile <- quantile(parse_obj@meta.data$percent.mt, seq(0,1,0.01))
    print(paste0("number of RNA count minimum: ", as.numeric(rna_quantile[2][[1]]), "; ", "number of RNA count maximum: ", as.numeric(rna_quantile[length(rna_quantile)-1][[1]])))
    print(paste0("number of maximum mt percentage: ", as.numeric(mt_quantile[length(mt_quantile)-1][[1]])))
    parse_obj <- subset(parse_obj, subset = percent.mt < as.numeric(mt_quantile[length(mt_quantile)-1][[1]]) & nCount_RNA > as.numeric(rna_quantile[2][[1]]) & nCount_RNA < as.numeric(rna_quantile[length(rna_quantile)-1][[1]]))


    return(parse_obj)
}



################################################################################################################
    # # name the transcripts without a gene name using gene_id
    # dataset <- "emory"
    # input <- "/scicomp/home-pure/rwo8/workdir/single_cell/emory/Processed/G24/" # do for each emory sample
    # output <- "/scicomp/home-pure/rwo8/workdir/single_cell/emory/work_RDS"
    # emory_genes <- read_tsv(file.path(input, "features.tsv"),col_names=FALSE)
    # colnames(emory_genes) <- c("gene_id", "gene_name", "gene_expression")
    # emory_genes_1 <- emory_genes %>% mutate(gene_name = mapply(function(x,y){
    #     ifelse(grepl("[.]", x), y, x)
    # }, gene_name, gene_id)) 
    # write_tsv(emory_genes_1,file.path(input, "features_1.tsv"),col_names=FALSE) 
    # # this new feature.tsv file need to be gzipped in their correspondinng folders 
################################################################################################################


# dataset <- "emory"
# input <- "/scicomp/home-pure/rwo8/workdir/single_cell/emory/Processed/Conts1/"
# output <- "/scicomp/home-pure/rwo8/workdir/single_cell/emory/work_RDS"
read_10x <- function(input){

    mtx_file <- Read10X(data.dir =input)
    # now create a seurat object with the data we read in
    emory_obj <- CreateSeuratObject(counts = mtx_file,project = basename(input),  min.genes = 1, min.cells = 1)

    emory_count <- emory_obj@assays$RNA@counts 
    emory_vir <- as.data.frame(t(as.data.frame(emory_count) %>% subset(rownames(.) %in% c("2012741656", "3000378993")))) %>% subset(rowSums(.) > 0) %>% rownames_to_column()
    emory_vir$sample <- basename(input)
    write.csv(emory_vir,file.path("/scicomp/home-pure/rwo8/workdir/single_cell/emory/output/tables", paste0(basename(input),"_virus_gene_expression_cells.csv" )))

    emory_obj@meta.data$orig.ident <- factor(sapply(emory_obj@meta.data$orig.ident, function(x){
    s <- switch(as.character(x),
            Conts1 = "NI",
            G23 = "GII3",
            G24 = "GII4"
            )
    paste(s, "_Emory", sep="_")
    }))

    emory_obj@meta.data$dataset <- factor(rep("Emory", nrow(emory_obj@meta.data)))

    Idents(emory_obj) <- emory_obj@meta.data$orig.ident

    emory_obj@meta.data$Sample <- factor(sapply(emory_obj@meta.data$orig.ident, function(x){
    unlist(strsplit(as.character(x), "_", fixed=TRUE))[1]
    }))

    # percentage of mitohondria reads per cell
    emory_obj[["percent.mt"]] <- PercentageFeatureSet(emory_obj , pattern = "^MT-")
    # percentage of ribosomomal reads per cell
    emory_obj[["percent.rb"]] <- PercentageFeatureSet(emory_obj, pattern = "^RP[SL]")

    rna_quantile <- quantile(emory_obj@meta.data$nCount_RNA, seq(0,1,0.01))
    mt_quantile <- quantile(emory_obj@meta.data$percent.mt, seq(0,1,0.01))
    print(paste0("number of RNA count minimum: ", as.numeric(rna_quantile[2][[1]]), "; ", "number of RNA count maximum: ", as.numeric(rna_quantile[length(rna_quantile)-1][[1]])))
    print(paste0("number of maximum mt percentage: ", as.numeric(mt_quantile[length(mt_quantile)-1][[1]])))
    emory_obj <- subset(emory_obj, subset = percent.mt < as.numeric(mt_quantile[length(mt_quantile)-1][[1]]) & nCount_RNA > as.numeric(rna_quantile[2][[1]]) & nCount_RNA < as.numeric(rna_quantile[length(rna_quantile)-1][[1]]))

    return(emory_obj)
}

if(dataset == "parse"){
    print("reading the parse sample into a seurat object ...")
    seurat_obj <- read_parse(input)
}else if(dataset == "emory"){
    print("reading the parse sample into a seurat object ...")
    seurat_obj <- read_10x(input)
}


###########################################################################################################
#
# normalize and run all seurat analysis
#
###########################################################################################################

print("normalizing and other seurat analysis...")
STC <- SCTransform(seurat_obj, method = "glmGamPoi",  vst.flavor = "v2", vars.to.regress = c("percent.mt"))

STC <- STC %>%
  RunPCA(verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE,metric = 'correlation') %>%
  RunTSNE(dims=1:20,dim.embed = 3) %>%
    FindNeighbors(reduction = "pca", dims =1:20, verbose = FALSE) %>%
    FindClusters(verbose = FALSE) # 3k cells normally has a good cluster with res 0.4 -1.2, higher res for bigger dataset

###########################################################################################################
#
# Find doublet
#
###########################################################################################################

print("find doublet in the current sample...")
library(SeuratDisk, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")
library(DoubletFinder, lib.loc="/scicomp/home-pure/rwo8/miniconda3/envs/python3/lib/R/library")


## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.STC <- paramSweep_v3(STC, PCs =1:20, sct = TRUE)
sweep.stats.STC <- summarizeSweep(sweep.res.STC, GT = FALSE)
bcmvn_STC <- find.pK(sweep.stats.STC)
pk <- as.numeric(as.character(bcmvn_STC[bcmvn_STC$BCmetric == max(bcmvn_STC$BCmetric),"pK"]))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(STC@meta.data$seurat_clusters)
nExp_poi <- round(0.03*nrow(STC@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
STC <- doubletFinder_v3(STC, PCs = 1:20, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

SaveH5Seurat(STC, output,overwrite = TRUE)



