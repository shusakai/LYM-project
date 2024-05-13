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


### Load the data ------------------------------
# The preprocessed data from subclustering.R
sc_T_cells <- readRDS("/path/to/sc_T_cells_hamaya.rds")


### CD8 pre vs during -------------------------------------------
sc_T_cells <- SetIdent(sc_T_cells, value = "celltype_pred_int_sub")
sc_CD8_cells <- subset(x = sc_T_cells, idents = c("Tregs", "CD4+ Tcm"), invert = TRUE)
sc_CD8_cells$RT.tag <- factor(sc_CD8_cells$RT.tag, levels = c("Pre","During","JustFinish","After"))

# DEG analysis, pre vs. During
sc_CD8_cells@active.assay <- "RNA"
Idents(sc_CD8_cells) <- "RT.tag"

de_markers.raw <- FindMarkers(sc_CD8_cells, 
                              ident.1 = "During", 
                              ident.2 = "Pre", 
                              assay = "RNA")

de_markers.raw$Significant <- ifelse(de_markers.raw$p_val_adj < 0.001 & abs(de_markers.raw$avg_log2FC) >= 0.5, 
                                     ifelse(de_markers.raw$avg_log2FC > 0.5, "Enriched in During", "Enriched in Pre"), "Intermediate")
de_markers.raw$Gene <- rownames(de_markers.raw)

# Volcano plot of DEG
# gene list for volcanoplot
d <- read.table('/path/to/volcano_label_gene.txt', header = F)
de_markers.raw.gene <- de_markers.raw %>% dplyr::filter(., Gene %in% d$V1)

# volcanoplot
voc2 <- ggplot(
  de_markers.raw, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant), size = 2) +
  scale_color_manual(values = c("firebrick","navy","grey")) +
  # annotation
  geom_text_repel(
    data = subset(de_markers.raw.gene, p_val_adj < 0.001 & abs(de_markers.raw.gene$avg_log2FC) >= 0.5),
    aes(label = Gene), max.overlaps = Inf, 
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  # hline
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
  # axis
  labs(x="log2(fold change)",
       y="-log10(p-value)") +
  # legend
  theme_classic() +
  theme(legend.position = "bottom") + 
  coord_cartesian(xlim = c(-1, 1), ylim = c(-10, 170))
ggsave(plot = voc2, filename = "/path/to/CD8_pre_during_volcano.png", height = 5, width = 5)

# up-DEG enrichment (Enriched in During)
up.deg <- subset(de_markers.raw, p_val_adj < 0.001 & de_markers.raw$avg_log2FC >= 0.5, select = Gene) 
up.deg <- bitr(up.deg$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

up.ego.bp <- enrichGO(gene = up.deg$ENTREZID, OrgDb = org.Hs.eg.db,
                      ont= "BP",pAdjustMethod = "BH",
                      pvalueCutoff= 0.05,
                      qvalueCutoff= 0.2,
                      readable= TRUE)

# Remove reduntant
up.ego.bp2 <- simplify(up.ego.bp)
bar1 <- barplot(up.ego.bp2, showCategory = 14,title="Enriched in During", x = "GeneRatio")
ggsave(plot = bar1, filename = "/path/to/CD8_during_pathway.png", height = 6, width = 12)


### CD8 pre vs just after -------------------------------------------
de_markers.raw <- FindMarkers(sc_CD8_cells, 
                              ident.1 = "JustFinish", 
                              ident.2 = "Pre", 
                              assay = "RNA")

de_markers.raw$Significant <- ifelse(de_markers.raw$p_val_adj < 0.001 & abs(de_markers.raw$avg_log2FC) >= 0.5, 
                                     ifelse(de_markers.raw$avg_log2FC > 0.5, "Enriched in Just After", "Enriched in Pre"), "Intermediate")
de_markers.raw$Gene <- rownames(de_markers.raw)

## Volcano plot of DEG
d <- read.table('/path/to/volcano_label_gene.txt', header = F)
de_markers.raw.gene <- de_markers.raw %>% dplyr::filter(., Gene %in% d$V1)

voc2 <- ggplot(
  de_markers.raw, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant), size = 2) +
  scale_color_manual(values = c("firebrick","navy","grey")) +
  # annotation
  geom_text_repel(
    data = subset(de_markers.raw.gene, p_val_adj < 0.05 & abs(de_markers.raw.gene$avg_log2FC) >= 0.5),
    aes(label = Gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  # hline
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
  # axis
  labs(x="log2(fold change)",
       y="-log10(p-value)") +
  # legend
  theme_classic() +
  theme(legend.position = "bottom") + 
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-10, 400))
ggsave(plot = voc2, filename = "/path/to/CD8_pre_justafter_volcano.png", height = 5, width = 5)

# up-DEG enrichment (Enriched in During)
up.deg <- subset(de_markers.raw, p_val_adj < 0.001 & de_markers.raw$avg_log2FC >= 0.5, select = Gene) 
up.deg <- bitr(up.deg$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

up.ego.bp <- enrichGO(gene = up.deg$ENTREZID, OrgDb = org.Hs.eg.db,
                      ont= "BP",pAdjustMethod = "BH",
                      pvalueCutoff= 0.05,
                      qvalueCutoff= 0.2,
                      readable= TRUE)

# Remove reduntant
up.ego.bp2 <- simplify(up.ego.bp)
bar1 <- barplot(up.ego.bp2, showCategory = 14,title="Enriched in Just after", x = "GeneRatio")
ggsave(plot = bar1, filename = "/path/to/CD8_justafter_pathway.png", height = 6, width = 12)


### CD8 pre vs after -------------------------------------------
de_markers.raw <- FindMarkers(sc_CD8_cells, 
                              ident.1 = "After", 
                              ident.2 = "Pre", 
                              assay = "RNA")

de_markers.raw$Significant <- ifelse(de_markers.raw$p_val_adj < 0.001 & abs(de_markers.raw$avg_log2FC) >= 0.5, 
                                     ifelse(de_markers.raw$avg_log2FC > 0.5, "Enriched in After", "Enriched in Pre"), "Intermediate")
de_markers.raw$Gene <- rownames(de_markers.raw)

## Volcano plot of DEG
d <- read.table('/path/to/volcano_label_gene.txt', header = F)
de_markers.raw.gene <- de_markers.raw %>% dplyr::filter(., Gene %in% d$V1)

voc2 <- ggplot(
  de_markers.raw, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant), size = 2) +
  scale_color_manual(values = c("firebrick","navy","grey")) +
  # annotation
  geom_text_repel(
    data = subset(de_markers.raw.gene, p_val_adj < 0.05 & abs(de_markers.raw.gene$avg_log2FC) >= 0.5),
    aes(label = Gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  # hline
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
  # axis
  labs(x="log2(fold change)",
       y="-log10(p-value)") +
  # legend
  theme_classic() +
  theme(legend.position = "bottom") + 
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-10, 400))
ggsave(plot = voc2, filename = "/path/to/CD8_pre_after_volcano.png", height = 5, width = 5)

# up-DEG enrichment (Enriched in During)
up.deg <- subset(de_markers.raw, p_val_adj < 0.001 & de_markers.raw$avg_log2FC >= 0.5, select = Gene) 
up.deg <- bitr(up.deg$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

up.ego.bp <- enrichGO(gene = up.deg$ENTREZID, OrgDb = org.Hs.eg.db,
                      ont= "BP",pAdjustMethod = "BH",
                      pvalueCutoff= 0.05,
                      qvalueCutoff= 0.2,
                      readable= TRUE)

# Remove reduntant
up.ego.bp2 <- simplify(up.ego.bp)
bar1 <- barplot(up.ego.bp2, showCategory = 14,title="Enriched in After", x = "GeneRatio")
ggsave(plot = bar1, filename = "/path/to/CD8_after_pathway.png", height = 6, width = 12)


### CD4 pre vs during -------------------------------------------
sc_T_cells <- SetIdent(sc_T_cells, value = "celltype_pred_int_sub")
sc_CD4_cells <- subset(x = sc_T_cells, idents = c("Tregs", "CD8+ Tcm"), invert = TRUE)
sc_CD4_cells$RT.tag <- factor(sc_CD4_cells$RT.tag, levels = c("Pre","During","JustFinish","After"))

# DEG analysis, pre vs. During
sc_CD4_cells@active.assay <- "RNA"
Idents(sc_CD4_cells) <- "RT.tag"

de_markers.raw <- FindMarkers(sc_CD4_cells, 
                              ident.1 = "During", 
                              ident.2 = "Pre", 
                              assay = "RNA")

de_markers.raw$Significant <- ifelse(de_markers.raw$p_val_adj < 0.001 & abs(de_markers.raw$avg_log2FC) >= 0.5, 
                                     ifelse(de_markers.raw$avg_log2FC > 0.5, "Enriched in During", "Enriched in Pre"), "Intermediate")
de_markers.raw$Gene <- rownames(de_markers.raw)

# Volcano plot of DEG
# gene list for volcanoplot
d <- read.table('/path/to/volcano_label_gene.txt', header = F)
de_markers.raw.gene <- de_markers.raw %>% dplyr::filter(., Gene %in% d$V1)

# volcanoplot
voc2 <- ggplot(
  de_markers.raw, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant), size = 2) +
  scale_color_manual(values = c("firebrick","navy","grey")) +
  # annotation
  geom_text_repel(
    data = subset(de_markers.raw.gene, p_val_adj < 0.001 & abs(de_markers.raw.gene$avg_log2FC) >= 0.5),
    aes(label = Gene), max.overlaps = Inf, 
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  # hline
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
  # axis
  labs(x="log2(fold change)",
       y="-log10(p-value)") +
  # legend
  theme_classic() +
  theme(legend.position = "bottom") + 
  coord_cartesian(xlim = c(-1, 1), ylim = c(-10, 170))
ggsave(plot = voc2, filename = "/path/to/CD4_pre_during_volcano.png", height = 5, width = 5)

# up-DEG enrichment (Enriched in During)
up.deg <- subset(de_markers.raw, p_val_adj < 0.001 & de_markers.raw$avg_log2FC >= 0.5, select = Gene) 
up.deg <- bitr(up.deg$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

up.ego.bp <- enrichGO(gene = up.deg$ENTREZID, OrgDb = org.Hs.eg.db,
                      ont= "BP",pAdjustMethod = "BH",
                      pvalueCutoff= 0.05,
                      qvalueCutoff= 0.2,
                      readable= TRUE)

# Remove reduntant
up.ego.bp2 <- simplify(up.ego.bp)
bar1 <- barplot(up.ego.bp2, showCategory = 14,title="Enriched in During", x = "GeneRatio")
ggsave(plot = bar1, filename = "/path/to/CD4_during_pathway.png", height = 6, width = 12)


### CD4 pre vs just after -------------------------------------------
de_markers.raw <- FindMarkers(sc_CD4_cells, 
                              ident.1 = "JustFinish", 
                              ident.2 = "Pre", 
                              assay = "RNA")

de_markers.raw$Significant <- ifelse(de_markers.raw$p_val_adj < 0.001 & abs(de_markers.raw$avg_log2FC) >= 0.5, 
                                     ifelse(de_markers.raw$avg_log2FC > 0.5, "Enriched in Just After", "Enriched in Pre"), "Intermediate")
de_markers.raw$Gene <- rownames(de_markers.raw)

## Volcano plot of DEG
d <- read.table('/path/to/volcano_label_gene.txt', header = F)
de_markers.raw.gene <- de_markers.raw %>% dplyr::filter(., Gene %in% d$V1)

voc2 <- ggplot(
  de_markers.raw, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant), size = 2) +
  scale_color_manual(values = c("firebrick","navy","grey")) +
  # annotation
  geom_text_repel(
    data = subset(de_markers.raw.gene, p_val_adj < 0.05 & abs(de_markers.raw.gene$avg_log2FC) >= 0.5),
    aes(label = Gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  # hline
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
  # axis
  labs(x="log2(fold change)",
       y="-log10(p-value)") +
  # legend
  theme_classic() +
  theme(legend.position = "bottom") + 
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-10, 400))
ggsave(plot = voc2, filename = "/path/to/CD4_pre_justafter_volcano.png", height = 5, width = 5)

# up-DEG enrichment (Enriched in During)
up.deg <- subset(de_markers.raw, p_val_adj < 0.001 & de_markers.raw$avg_log2FC >= 0.5, select = Gene) 
up.deg <- bitr(up.deg$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

up.ego.bp <- enrichGO(gene = up.deg$ENTREZID, OrgDb = org.Hs.eg.db,
                      ont= "BP",pAdjustMethod = "BH",
                      pvalueCutoff= 0.05,
                      qvalueCutoff= 0.2,
                      readable= TRUE)

# Remove reduntant
up.ego.bp2 <- simplify(up.ego.bp)
bar1 <- barplot(up.ego.bp2, showCategory = 14,title="Enriched in Just after", x = "GeneRatio")
ggsave(plot = bar1, filename = "/path/to/CD4_justafter_pathway.png", height = 6, width = 12)


### CD4 pre vs after -------------------------------------------
de_markers.raw <- FindMarkers(sc_CD4_cells, 
                              ident.1 = "After", 
                              ident.2 = "Pre", 
                              assay = "RNA")

de_markers.raw$Significant <- ifelse(de_markers.raw$p_val_adj < 0.001 & abs(de_markers.raw$avg_log2FC) >= 0.5, 
                                     ifelse(de_markers.raw$avg_log2FC > 0.5, "Enriched in After", "Enriched in Pre"), "Intermediate")
de_markers.raw$Gene <- rownames(de_markers.raw)

## Volcano plot of DEG
d <- read.table('/path/to/volcano_label_gene.txt', header = F)
de_markers.raw.gene <- de_markers.raw %>% dplyr::filter(., Gene %in% d$V1)

voc2 <- ggplot(
  de_markers.raw, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significant), size = 2) +
  scale_color_manual(values = c("firebrick","navy","grey")) +
  # annotation
  geom_text_repel(
    data = subset(de_markers.raw.gene, p_val_adj < 0.05 & abs(de_markers.raw.gene$avg_log2FC) >= 0.5),
    aes(label = Gene),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  # hline
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1,lty=4,col="black",lwd=0.8) +
  # axis
  labs(x="log2(fold change)",
       y="-log10(p-value)") +
  # legend
  theme_classic() +
  theme(legend.position = "bottom") + 
  coord_cartesian(xlim = c(-1.2, 1.2), ylim = c(-10, 400))
ggsave(plot = voc2, filename = "/path/to/CD4_pre_after_volcano.png", height = 5, width = 5)

# up-DEG enrichment (Enriched in During)
up.deg <- subset(de_markers.raw, p_val_adj < 0.001 & de_markers.raw$avg_log2FC >= 0.5, select = Gene) 
up.deg <- bitr(up.deg$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

up.ego.bp <- enrichGO(gene = up.deg$ENTREZID, OrgDb = org.Hs.eg.db,
                      ont= "BP",pAdjustMethod = "BH",
                      pvalueCutoff= 0.05,
                      qvalueCutoff= 0.2,
                      readable= TRUE)

# Remove reduntant
up.ego.bp2 <- simplify(up.ego.bp)
bar1 <- barplot(up.ego.bp2, showCategory = 14,title="Enriched in After", x = "GeneRatio")
ggsave(plot = bar1, filename = "/path/to/CD4_after_pathway.png", height = 6, width = 12)

