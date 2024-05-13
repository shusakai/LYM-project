### Setup the Seurat object --------------------------------------------
# Load Packages
packages <- c("Seurat","glmGamPoi","DOSE","ggplot2","ggrepel","ggpubr","grid","cowplot","patchwork",
              "DESeq2","wesanderson","pheatmap","dplyr",
              "SingleR", "celldex",
              "org.Hs.eg.db","clusterProfiler", "escape","dittoSeq",
              "ComplexHeatmap","FlexDotPlot","tidyverse")
installed_packages <- packages %in% rownames(installed.packages())  # Install packages not yet installed
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
}

invisible(lapply(packages, library, character.only = TRUE)) # Packages loading

# load hamaya-san data
sc_T_cells <- readRDS("path/to/sc_T_cells_hamaya.rds")
sc_Epithelial <- readRDS("path/to/sc_Epithelial_hamaya.rds")

# load du-san data
sc_major <- readRDS("/Volumes/Carrot/home/data/scRNA/sc_major.rds")
sc_Myeloid <- subset(x = sc_major, idents = "Myeloid_cells")
sc_Myeloid$RT.tag <- factor(sc_Myeloid$RT.tag, levels = c("Pre","During","JustFinish","After"))

sc_Myeloid@active.assay <- "integrated"
sc_Myeloid <- RunPCA(sc_Myeloid, features = VariableFeatures(object = sc_Myeloid))
ElbowPlot(sc_Myeloid)
PC.num = 1:15

sc_Myeloid <- RunUMAP(sc_Myeloid, dims = PC.num)
sc_Myeloid<- FindNeighbors(sc_Myeloid, dims = PC.num)
sc_Myeloid <- FindClusters(sc_Myeloid, resolution = 0.5)
SingleR_ref_hpca <- HumanPrimaryCellAtlasData(ensembl = FALSE)
SingleR_ref_bpe <- BlueprintEncodeData(ensembl = FALSE)
testdata <- GetAssayData(sc_Myeloid,slot = "data")
clusters <- sc_Myeloid@meta.data$seurat_clusters
cellpred <- SingleR(testdata,
                    ref = SingleR_ref_bpe, 
                    labels = SingleR_ref_bpe$label.fine,
                    clusters = clusters, 
                    method = "cluster",
                    assay.type.test = "logcounts",
                    assay.type.ref = "logcounts")
celltype <- data.frame(ClusterID = rownames(cellpred), celltype = cellpred$labels, stringsAsFactors = F)
sc_Myeloid@meta.data$celltype_pred_int_sub <- "NA"
for (id in 1:nrow(celltype)) {
  sc_Myeloid@meta.data[which(sc_Myeloid@meta.data$seurat_clusters == celltype$ClusterID[id]),"celltype_pred_int_sub"] <- celltype$celltype[id]
}

unique(sc_Myeloid$celltype_pred_int_sub)
Idents(sc_Myeloid) <- sc_Myeloid$celltype_pred_int_sub

reduct<- "umap" # or "tsne"
p0 <- DimPlot(sc_Myeloid, 
              group.by = "seurat_clusters",
              reduction = reduct,
              repel = TRUE,
              label = T,
              shuffle = T) + ggtitle("Myeloid cells subcluster") + NoLegend()

p1 <- DimPlot(sc_Myeloid, 
              group.by = "celltype_pred_int_sub",
              order = c("DC","Monocytes","Macrophages M1","Macrophages","Neutrophils"),
              reduction = reduct,
              repel = TRUE,
              label = T,
              shuffle = T) + ggtitle("Myeloid cells subtype") + NoLegend()

# save output
plot_grid(p0, p1, ncol=2, align="hv")
Idents(sc_Myeloid) <- "celltype_pred_int_sub"
sc_Myeloid <- RenameIdents(object = sc_Myeloid, `Macrophages M1` = "Macrophages")
sc_Myeloid$celltype_pred_int_sub <- Idents(sc_Myeloid)
# save
#saveRDS(sc_Myeloid, "/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/sc_Myeloid_cells_du.rds")

# read Myeloid cells
#sc_Myeloid <- readRDS("/Volumes/Project/home/sakai/20230704_Tcell_Epithelial_scRNA_project/data/sc_Myeloid_cells_du.rds")

# rename indents
sc_T_cells <- RenameIdents(object = sc_T_cells, `CD8+ Tem` = "CD8+ T cells", `CD4+ Tcm` = "CD4+ T cells")
sc_T_cells$celltype_pred_int_sub <- Idents(sc_T_cells)

sc_Myeloid@active.assay <- "RNA"

# merge two object
sc_T_Epithelial_cells <- merge(sc_T_cells, y=sc_Epithelial)
sc_T_Epithelials_Myeloid_cells <- merge(sc_Myeloid, y=sc_T_Epithelial_cells)

# factorize
sc_T_Epithelials_Myeloid_cells$RT.tag <- factor(sc_T_Epithelials_Myeloid_cells$RT.tag, levels = c("Pre", "During", "JustFinish","After"))

### Draw dotplot of interested gene sets(Fig.4B) ####
many_dotplots <- function(gene_signatures, file_name, celltype){
  one_dotplot <- function(signature,celltype){
    
    # Prepare the dataset using Seurat
    dot1 <- DotPlot(celltype,
                    split.by = "RT.tag",
                    assay = "RNA",
                    group.by = "celltype_man",
                    features = signature,
                    cols = 1:8,
                    dot.scale = 10) + xlab("") + ylab("")
    
    df_dotplot <- dot1$data %>% dplyr::select(-colors) %>% filter(pct.exp > 0)
    rownames(df_dotplot) <- 1:nrow(df_dotplot)
    
    # Separate cell and time tag
    df_dotplot$cell_time <- df_dotplot$id
    df_dotplot <- tidyr::separate(data = df_dotplot, col = id, into = c("cell.tag","time.tag"),sep = "_") %>% 
      dplyr::select(features.plot,cell_time,everything())
    
    # Order time tag
    df_dotplot$time.tag <- factor(df_dotplot$time.tag, levels = c("Pre", "During", "JustFinish","After"))
    
    df_dotplot$cell.tag <- factor(df_dotplot$cell.tag, levels = c("Macrophages", "Neutrophils", "Monocytes", "DC", "CD8+ T cells", "CD4+ T cells", "Tregs","Epithelial cells"))
    
    # Draw dotplot
    dot2 <- ggplot(data = df_dotplot) + 
      geom_point(mapping = aes(x = time.tag, 
                               y = features.plot, 
                               size = pct.exp, 
                               color = avg.exp.scaled)) +
      scale_colour_gradient2() + 
      scale_size(range = c(0,6),
                 name = "%Expression",
                 breaks = c(1,20,40,60,80,100)) +
      
      # add facet by cell type
      facet_grid(cols = vars(cell.tag)) + theme_bw() + NoGrid() +
      theme(axis.title = element_blank(),
            axis.text.y = element_text(size = 10,color = "black"),
            axis.text.x = element_text(size = 8,angle = 45,hjust = 1,colour = "black"),
            legend.title = element_text(size = 8),
            legend.text =  element_text(size = 7),
            legend.key.width=unit(0.3,"cm"),
            legend.key.height=unit(0.3,"cm"),
            legend.spacing = unit(0.1,"cm"))
    
  }
  for (signature in names(gene_signatures)) {
    # change signature as a vector(Seurat parameter)
    sig <- get(signature) 
    
    # Draw one dotplot
    p1 <- one_dotplot(sig,celltype)
    
    # Save the dotplot
    new_file_name <- paste0(file_name, "5.2.1_dot_",signature, ".png")
    print(new_file_name)
    height = length(levels(p1$data$features.plot))
    weight = length(levels(p1$data$cell_time))
    ggsave(p1, filename = new_file_name, width = 0.45*weight, height = 0.35*height)
    
    # reset the file name
    new_file_name <- NULL
  }
}

# IFN asocciated gene
IFNG <- c(
  "IRF1", "IRF9", # downstream gene
  "OAS2", "MX1", 
  "IFNGR1", "IFNGR2", 
  "IFNG" # IFNG
)

# 20230808 updated: immune checkpoint
ICI <- c(
  "LAG3", "LGALS3", "CLEC4G", # LAG3
  "BTLA", "TNFRSF14", # BTLA
  "CD244","CD48", # CD244
  "PVR", "TIGIT", # TIGIT
  "CD80", "CD86", "CTLA4",  # CD80
  "PDCD1", "PDCD1LG2", "CD274" # PD-1
)

# 20230808 updated: TGFB
TGFB <- rev(c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3", "SMAD2", "SMAD3", "CDKN1A", "MMP9", "MYC"))


# Full gene set
gene_signatures = list(IFNG = IFNG,
                       Chemokines = Chemokines,
                       ICI = ICI, 
                       ligand_receptor = ligand_receptor, 
                       TGFB=TGFB, 
                       TEST=TEST, 
                       LTR=LTR)

subset <- "many_dot" 
output_dir <- "path/to/output"

rownames(sc_T_Epithelials_Myeloid_cells@assays$RNA)

sc_T_Epithelials_Myeloid_cells$celltype_man <- sc_T_Epithelials_Myeloid_cells$celltype_pred_int_sub
sc_T_Epithelials_Myeloid_cells$celltype_man <- factor(sc_T_Epithelials_Myeloid_cells$celltype_man, levels = c(
  "Macrophages", "Neutrophils", "Monocytes", "DC", "CD8+ T cells", "CD4+ T cells", "Tregs","Epithelial cells"
))
Idents(sc_T_Epithelials_Myeloid_cells) <- sc_T_Epithelials_Myeloid_cells$celltype_man


many_dotplots(gene_signatures, file_name = paste0(output_dir,"/",subset,"/"), sc_T_Epithelials_Myeloid_cells)
# save
saveRDS(sc_T_Epithelials_Myeloid_cells, "path/to/output/sc_Mye_T_Epi.rds")


### separate time cource of seurat object -----------------------------------------------
unique(Idents(sc_T_Epithelials_Myeloid_cells))
unique(sc_T_Epithelials_Myeloid_cells$celltype_pred_int_sub)

# pre
sc_T_Epithelials_Myeloid_cells_pre <- subset(
  sc_T_Epithelials_Myeloid_cells, 
  idents = c("Pre")
)
# during
sc_T_Epithelials_Myeloid_cells_during <- subset(
  sc_T_Epithelials_Myeloid_cells, 
  idents = c("During")
)
# JustFinish
sc_T_Epithelials_Myeloid_cells_just <- subset(
  sc_T_Epithelials_Myeloid_cells, 
  idents = c("JustFinish")
)
# After
sc_T_Epithelials_Myeloid_cells_after <- subset(
  sc_T_Epithelials_Myeloid_cells, 
  idents = c("After")
)

# save
saveRDS(sc_T_Epithelials_Myeloid_cells_pre, 
        "path/to/sc_pre.rds")
saveRDS(sc_T_Epithelials_Myeloid_cells_during, 
        "path/to/sc_during.rds")
saveRDS(sc_T_Epithelials_Myeloid_cells_just, 
        "path/to/sc_justfinish.rds")
saveRDS(sc_T_Epithelials_Myeloid_cells_after, 
        "path/to/sc_after.rds")