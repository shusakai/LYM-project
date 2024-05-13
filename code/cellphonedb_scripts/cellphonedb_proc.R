# load packages -----------------------------------------------------------
require(tibble)
require(biomaRt)
require(tidyr)
require(dplyr)
require(Seurat)

# define function -----------------------------------------------------------
shaping_fun <- function(obj, name_, out_path){
  # shaping
  alldata <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  alldata@active.assay <- "RNA"
  allgenes <- rownames(alldata)
  matrix1 <- as.data.frame(alldata@assays$RNA@data)
  matrix1 <- matrix1[rowSums(matrix1[,2:dim(matrix1)[2]])!=0,]
  
  ### If you are using a mouse data, then its needed to convert the gene names to human orthologs
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")
  genesV2 <- getBM(attributes= c("ensembl_gene_id","hgnc_symbol"), filters= "hgnc_symbol", values = rownames(alldata@assays$RNA@data), mart= human, uniqueRows=T)
  #matrix1 <- matrix1[match(genesV2$hgnc_symbol, rownames(alldata), nomatch=F),]
  #matrix1$gene <- genesV2[match(genesV2$hgnc_symbol, rownames(alldata), nomatch=F),]$ensembl_gene_id
  matrix1$hgnc_symbol <- rownames(matrix1)
  matrix1_joined <- merge(matrix1, genesV2, by = "hgnc_symbol", all.x = T, all.y = F)
  matrix1_joined <- matrix1_joined[is.na(matrix1_joined$ensembl_gene_id) == F, ]
  matrix1_joined <- matrix1_joined %>% distinct(hgnc_symbol, .keep_all = T) 
  rownames(matrix1_joined) <- matrix1_joined$hgnc_symbol
  matrix1_joined <- matrix1_joined[, colnames(matrix1_joined) != "hgnc_symbol"]
  
  ## If the cluster names are categorical, you will need to convert it to numerical
  alldata@meta.data$celltype_pred_int_sub <- as.factor(alldata@meta.data$celltype_pred_int_sub)
  metadata <- data.frame(cells=colnames(matrix1), cluster=alldata@meta.data[colnames(matrix1), "celltype_pred_int_sub"])
  
  # saves
  write.table(matrix1_joined, paste(out_path, paste(name_, "_filtered_hcount.txt", sep=""), sep="/"), row.names=T,sep = "\t")
  write.table(metadata, paste(out_path, paste(name_, "_filtered_meta.txt", sep=""), sep="/"), row.names=FALSE, sep = "\t")
}

# main -----------------------------------------------------------
# obtain argument
args <- commandArgs(trailingOnly = T)
obj_path <- args[1]
proj_name <- args[2]
out_path <- args[3]

# load seurat object
seu_obj <- readRDS(obj_path)

# run 
shaping_fun(seu_obj, proj_name, out_path)

