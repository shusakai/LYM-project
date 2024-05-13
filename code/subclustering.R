### Setup the Seurat-------------------------------------------------------
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

# Load reference data set
SingleR_ref_hpca <- HumanPrimaryCellAtlasData(ensembl = FALSE)
SingleR_ref_bpe <- BlueprintEncodeData(ensembl = FALSE)

### Load the single cell RNA sequence count data --------------------------
# Oyoshi et al Sci. Adv. 2023
sc_major <- readRDS("/path/to/sc_major.rds")

### Subclustering of subclustering------------------------------
subset <- "T_cells" 
sc_T_cells <- subset(x = sc_major, idents = "T_cells")
sc_T_cells@active.assay <- "integrated"

png("/path/to/RT_time_c_T_cells.png", width = 400, height = 300)
barplot(table(sc_T_cells$RT.tag))
dev.off() 
reduct <-"umap"
#reduct <-"tsne"

# DimPlot
p0 <- DimPlot(sc_T_cells, 
              group.by = "seurat_clusters",label.size = 6, 
              reduction = reduct,
              label = T,
              shuffle = T) + ggtitle("T cell subclusters") + NoLegend()

p1 <- DimPlot(sc_T_cells, 
              group.by = "celltype_man2",label.size = 6, 
              label = T,
              shuffle = T) + ggtitle("T cell subset") + NoLegend()
RT_dim <- DimPlot(sc_T_cells, label.size = 6, 
                  group.by = "RT.tag",
                  reduction = reduct,
                  label = T, repel = TRUE,
                  shuffle = T) + ggtitle("T cell RT label")

plot_grid(p0, p1, RT_dim, ncol=2, nrow=2, align="hv")
#ggsave(filename = paste(dataset,subset,"5.2.1_subdim_integrated_preCluster_T_cell.png", sep="/"), width = 20, height = 15)
sc_T_cells$old_seurat_clusters <- sc_T_cells$seurat_clusters


# Re-run SingleR on T_cells
sc_T_cells <- RunPCA(sc_T_cells, features = VariableFeatures(object = sc_T_cells))
ElbowPlot(sc_T_cells)

PC.num = 1:10 

sc_T_cells <- RunUMAP(sc_T_cells, dims = PC.num)
sc_T_cells <- RunTSNE(sc_T_cells, dims = PC.num)

sc_T_cells<- FindNeighbors(sc_T_cells, dims = PC.num)
sc_T_cells <- FindClusters(sc_T_cells, resolution = 0.5)

testdata <- GetAssayData(sc_T_cells,slot = "data")
clusters <- sc_T_cells@meta.data$seurat_clusters
cellpred <- SingleR(testdata,
                    ref = SingleR_ref_bpe, 
                    labels = SingleR_ref_bpe$label.fine,
                    clusters = clusters, 
                    method = "cluster",
                    assay.type.test = "logcounts",
                    assay.type.ref = "logcounts")
celltype <- data.frame(ClusterID = rownames(cellpred), celltype = cellpred$labels, stringsAsFactors = F)
sc_T_cells@meta.data$celltype_pred_int_sub <- "NA"
for (id in 1:nrow(celltype)) {
  sc_T_cells@meta.data[which(sc_T_cells@meta.data$seurat_clusters == celltype$ClusterID[id]),"celltype_pred_int_sub"] <- celltype$celltype[id]
}

p0 <- DimPlot(sc_T_cells, 
              group.by = "seurat_clusters",label.size = 6, 
              reduction = reduct,
              label = T,
              shuffle = T) + ggtitle("T_cells cells subclusters") + NoLegend()

p1 <- DimPlot(sc_T_cells, 
              group.by = "celltype_pred_int_sub",label.size = 6, 
              label = T,
              shuffle = T) + ggtitle("T_cells cells subset") + NoLegend()
RT_dim <- DimPlot(sc_T_cells, label.size = 6, 
                  group.by = "RT.tag",
                  reduction = reduct,
                  label = T, repel = TRUE, 
                  shuffle = T) + ggtitle("T_cells cells RT label")

plot_grid(p0, p1, RT_dim, ncol=2, nrow=2, align="hv")
#ggsave(filename = paste(dataset,subset,"dimplot_integrated_T_cells_singleR.png", sep="/"), width = 20, height = 15)

#detele mesangial cell
Idents(sc_T_cells) <- "celltype_pred_int_sub"
sc_T_cells_no_mes <- subset(sc_T_cells, ident = c("Mesangial cells"), invert=T)
unique(Idents(sc_T_cells_no_mes))

# singleR plot 
p1 <- DimPlot(sc_T_cells_no_mes, 
              group.by = "celltype_pred_int_sub",label.size = 6, 
              label = T,
              shuffle = T) + ggtitle("SingleR") + NoLegend()
plot_grid(p1, ncol=1, nrow=1, align="hv")
ggsave(filename = "/path/to/T_cells_annotation_Fig1C.png", width = 5, height = 5, dpi = 500)

### Marker gene plot ------------------------------------------------------------------
sc_T_cells_no_mes@active.assay <- "RNA"
FP <- wes_palette("Zissou1", 5, type = "discrete")

# CD8
g1 <- FeaturePlot(sc_T_cells_no_mes, "CD8A", reduction = reduct,
                  pt.size = 1, min.cutoff = "q10", max.cutoff = "q90")+
  NoAxes()+
  theme(plot.title = element_text(size=20, face="bold"),
        axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"))
plot_grid(g1, ncol=1, nrow=1, align="hv")
ggsave(filename = "/path/to/20230704_Tcell_Epithelial_scRNA_project/data/T_cells_CD8.png", width = 6, height = 5, dpi = 500)

# CD4
g1 <- FeaturePlot(sc_T_cells_no_mes, "CD4", reduction = reduct,
                  pt.size = 1, min.cutoff = "q10", max.cutoff = "q90")+
  NoAxes()+
  theme(plot.title = element_text(size=20, face="bold"),
        axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"))
plot_grid(g1, ncol=1, nrow=1, align="hv")
ggsave(filename = "/path/to/20230704_Tcell_Epithelial_scRNA_project/data/T_cells_CD4.png", width = 6, height = 5, dpi = 500)

# FOXP3
g1 <- FeaturePlot(sc_T_cells_no_mes, "FOXP3", reduction = reduct,
                  pt.size = 1, min.cutoff = "q10", max.cutoff = "q90")+
  NoAxes()+
  theme(plot.title = element_text(size=20, face="bold"),
        axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"))
plot_grid(g1, ncol=1, nrow=1, align="hv")
ggsave(filename = "/path/to/20230704_Tcell_Epithelial_scRNA_project/data/T_cells_Treg.png", width = 6, height = 5, dpi = 500)

# CD25
g1 <- FeaturePlot(sc_T_cells_no_mes, "IL2RA", reduction = reduct,
                  pt.size = 1, min.cutoff = "q10", max.cutoff = "q90")+
  NoAxes()+
  theme(plot.title = element_text(size=20, face="bold"),
        axis.text=element_text(size=14,face="bold"),
        axis.title=element_text(size=16,face="bold"))
plot_grid(g1, ncol=1, nrow=1, align="hv")
ggsave(filename = "/path/to/20230704_Tcell_Epithelial_scRNA_project/data/T_cells_ CD25.png", width = 6, height = 5, dpi = 500)


### Save the object ----------------------------------------
# remove NK cells 
sc_T_cells <- SetIdent(sc_T_cells, value = "celltype_pred_int_sub")
sc_T_cells<- subset(x = sc_T_cells, idents = c("NK cells","Mesangial cells"), invert = TRUE)
sc_T_cells$RT.tag <- factor(sc_T_cells$RT.tag, levels = c("Pre","During","JustFinish","After"))
# save
saveRDS(sc_T_cells, "/path/to/20230704_Tcell_Epithelial_scRNA_project/data/sc_T_cells_hamaya.rds")