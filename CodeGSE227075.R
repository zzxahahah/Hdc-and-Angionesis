reticulate::use_python('F:/anaconda/envs/SCP_env/python.exe', required = TRUE)

library(Seurat)
library(SingleR)
library(SeuratWrappers)
library(monocle3)
library(reticulate)
library(Matrix)
library(ggplot2)
library(cowplot)
library(magrittr)
library(patchwork)
library(dplyr)
library(hdf5r)
library(scuttle)
library(slingshot)
library(sceasy)
library(RColorBrewer)
library(BiocParallel)
library(decontX)
library(harmony)
library(corrplot)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(CytoTRACE)
library(SCP)
library(CellChat)
library(igraph)

register(SnowParam(workers = 6, progressbar = TRUE))
setwd("C:/Users/zjf/Desktop/skeleton muscle")

dataday1_dir <- './GSE227075'
GSE227075 <- Read10X(data.dir = dataday1_dir, gene.column = 2)
GSE227075 <- CreateSeuratObject(counts = GSE227075, project = "GSE227075", min.cells = 3, min.features = 200)

#QC
GSE227075 <- RunCellQC(srt = GSE227075, qc_metrics = c("doublets","outlier","mito","ribo","ribo_mito_ratio"), 
                       mito_threshold = 5, db_rate = 0.075)
Srt.GSE227075 <- subset(GSE227075, CellQC == "Pass")
table(GSE227075$CellQC)
rm(GSE227075)
gc()

#standard pipeline
Srt.GSE227075 <- NormalizeData(Srt.GSE227075, verbose = F)
Srt.GSE227075 <- ScaleData(Srt.GSE227075, verbose = F)
Srt.GSE227075 <- FindVariableFeatures(Srt.GSE227075, selection.method = "vst", nfeatures = 2000, verbose = F)
Srt.GSE227075 <- RunPCA(Srt.GSE227075, npcs = 50, verbose = F)
Srt.GSE227075 <- RunUMAP(Srt.GSE227075, reduction = "pca", dims = 1:30, verbose = F)

#batch annotation
cell_names <- Idents(Srt.GSE227075)

Batchorder = rep("0",length(cell_names))
Batchorder[grepl("sham", colnames(Srt.GSE227075))] = "0"
Batchorder[grepl("day1", colnames(Srt.GSE227075))] = "1"
Batchorder[grepl("day3", colnames(Srt.GSE227075))] = "2"
Batchorder[grepl("day7", colnames(Srt.GSE227075))] = "3"
names(Batchorder) = colnames(Srt.GSE227075)
table(Batchorder)

Batchid = rep("sham",length(cell_names))
Batchid[grepl("sham", colnames(Srt.GSE227075))] = "Sham"
Batchid[grepl("day1", colnames(Srt.GSE227075))] = "Day1"
Batchid[grepl("day3", colnames(Srt.GSE227075))] = "Day3"
Batchid[grepl("day7", colnames(Srt.GSE227075))] = "Day7"
names(Batchid) = colnames(Srt.GSE227075)
table(Batchid)

Srt.GSE227075 <- AddMetaData(object = Srt.GSE227075, metadata = Batchorder, col.name = "Batchorder")
table(Srt.GSE227075@meta.data$Batchorder)

Srt.GSE227075 <- AddMetaData(object = Srt.GSE227075, metadata = Batchid, col.name = "Batchid")
table(Srt.GSE227075@meta.data$Batchid)

#decontX
Srt.GSE227075.counts <- Srt.GSE227075@assays$RNA@counts
decontX_results <- decontX(Srt.GSE227075.counts)
Srt.GSE227075$Contamination = decontX_results$contamination
Srt.GSE227075 = Srt.GSE227075[,Srt.GSE227075$Contamination < 0.2]
DimPlot(Srt.GSE227075,reduction = "umap")

#total harmony
Srt.GSE227075 <- Srt.GSE227075 %>% 
  RunHarmony("Batchid", plot_convergence = T)
harmony_embeddings <- Embeddings(Srt.GSE227075, 'harmony')
harmony_embeddings[1:5, 1:5]
DimPlot(object = Srt.GSE227075, reduction = 'umap', group.by = "Batchid", pt.size=0.1)

#cell cycling evaluation
cc.genes <- readRDS("./mouse_cell_cycle_genes.rds")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

Srt.GSE227075 <- SetIdent(Srt.GSE227075, value = 'Celltype')
Srt.GSE227075 <- CellCycleScoring(Srt.GSE227075, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(Srt.GSE227075[[]])
RidgePlot(Srt.GSE227075, 
          features = c("Pcna", "Top2a", "Mcm6", "Mki67"), 
          ncol = 2)
Srt.GSE227075$CC.Difference <- Srt.GSE227075$S.Score - Srt.GSE227075$G2M.Score
Srt.GSE227075 <- ScaleData(Srt.GSE227075, vars.to.regress = "CC.Difference", features = rownames(Srt.GSE227075))

#save
saveRDS(Srt.GSE227075,"srtGSE227075.rds")
Srt.GSE227075 <- readRDS("srtGSE227075.rds")

#cluster
ElbowPlot(object = Srt.GSE227075, ndims = 50)
use.pcs = 1:35
plan("multisession", workers = 4)
Srt.GSE227075 <- FindNeighbors(object = Srt.GSE227075, reduction = "harmony", dims= use.pcs)
Srt.GSE227075 <- FindClusters(object = Srt.GSE227075, reduction = "harmony", dims.use = use.pcs, resolution = seq(0.1,2,0.1), print.output = FALSE, save.SNN = TRUE)
Srt.GSE227075 <- RunUMAP(object = Srt.GSE227075, reduction = "harmony", dims = use.pcs)
sapply(grep("^Srt.GSE227075_snn_res",colnames(Srt.GSE227075@meta.data),value = TRUE), function(x) length(unique(Srt.GSE227075@meta.data[,x])))
p1 <- CellDimPlot(srt = Srt.GSE227075, group.by = "RNA_snn_res.2", reduction = "UMAP", theme_use = "theme_blank")
p1
ggsave("Totalclusterres0.6.pdf",p1,width = 6, height = 4)
DimPlot(object = Srt.GSE227075, reduction = 'umap', group.by = "RNA_snn_res.2", pt.size=0.1, label = TRUE, repel = TRUE) + NoLegend()
Srt.GSE227075 <- SetIdent(Srt.GSE227075, value = 'Celltype')

plan("sequential")
Srt.GSE227075 <- NormalizeData(Srt.GSE227075, verbose = FALSE)
all.genes <- rownames(Srt.GSE227075)
Srt.GSE227075 <- ScaleData(Srt.GSE227075, features = all.genes, verbose = FALSE)
plan("multisession", workers = 6)
markers_all_RNA <- FindAllMarkers(object = Srt.GSE227075, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "LR")
write.csv(markers_all_RNA, file = "cluster_allmarker_Celltype.csv")
table(Srt.GSE227075@meta.data$RNA_snn_res.0.6,Srt.GSE227075@meta.data$Batchid)

#SingleR
mrsd.se <- MouseRNAseqData()
Srt.GSE227075.se <- as.SingleCellExperiment(Srt.GSE227075, assay = "RNA")
mrsd.common <- intersect(rownames(Srt.GSE227075.se), rownames(mrsd.se))
mrsd.se <- mrsd.se[mrsd.common,]
Srt.GSE227075.se <- Srt.GSE227075.se[mrsd.common,]
Srt.GSE227075.se <- logNormCounts(Srt.GSE227075.se)
clusters <- Srt.GSE227075@meta.data$RNA_snn_res.1.2 #change the resolution
cluster.sr <- SingleR(test = Srt.GSE227075.se, ref = mrsd.se, method = "cluster",clusters = clusters, labels = mrsd.se$label.main)
Celltype = data.frame(ClusterID=rownames(cluster.sr), Celltype=cluster.sr$labels, stringsAsFactors = FALSE)
head(Celltype)
write.csv(Celltype,"Celltype_singleR.csv",row.names = FALSE)

#original annotation
annotation.origin <- read.csv("originalannotation.csv",head = TRUE)
cell_names <- colnames(Srt.GSE227075)
Celltype = rep("0",length(cell_names))
annotation.origin <- subset(annotation.origin, cellID %in% cell_names)
Celltype <- annotation.origin$annotation[cell_names %in% annotation.origin$cellID]
names(Celltype) = colnames(Srt.GSE227075)
table(Celltype)
Srt.GSE227075 <- AddMetaData(object = Srt.GSE227075, metadata = Celltype, col.name = "Celltype")
table(Srt.GSE227075@meta.data$Celltype)
p2 <- CellDimPlot(srt = Srt.GSE227075, group.by = "Celltype", reduction = "UMAP", label = TRUE, theme_use = "theme_blank")
p2
ggsave("Celltypeplot.pdf",p2,width = 6, height = 4)

#annotation after manual check
table(Srt.GSE227075$RNA_snn_res.2)
Srt.GSE227075 <- SetIdent(Srt.GSE227075, value = 'RNA_snn_res.2')
Srt.GSE227075$Celltype <- NULL
new.cluster.ids <- c("FAPs",
                     "Macrophages",
                     "Macrophages",
                     "Macrophages",
                     "Macrophages",
                     "Macrophages",
                     "B cells",
                     "MuSCs",
                     "FAPs",
                     "Macrophages",
                     "Macrophages",
                     "FAPs",
                     "Neutrophils",
                     "Cycling basal cells",
                     "FAPs",
                     "FAPs",
                     "Macrophages",
                     "Macrophages",
                     "Macrophages",
                     "Enodothelial cells",
                     "MuSCs",
                     "FAPs",
                     "FAPs",
                     "Schwann cells",
                     "NK/T cells",
                     "B cells",
                     "Tenocytes",
                     "Enodothelial cells",
                     "FAPs",
                     "MuSCs",
                     "Macrophages",
                     "NK/T cells",
                     "Erythrocytes",
                     "Epithelial cells",
                     "FAPs",
                     "NK/T cells",
                     "MuSCs",
                     "Enodothelial cells",
                     "FAPs",
                     "Macrophages",
                     "B cells",
                     "FAPs",
                     "Erythrocytes",
                     "Macrophages",
                     "Neutrophils")
names(new.cluster.ids) <- levels(Srt.GSE227075) 
Srt.GSE227075 <- RenameIdents(Srt.GSE227075, new.cluster.ids)
Srt.GSE227075[["Celltype"]] <- Idents(Srt.GSE227075)
p2 <- CellDimPlot(srt = Srt.GSE227075, group.by = "Celltype", reduction = "UMAP", theme_use = "theme_blank")
p2
ggsave("Celltypeplot.pdf",p2,width = 6, height = 4)


macrophage <- FeatureDimPlot(Srt.GSE227075, features = c("Arg1","Msr1","Mrc1"), theme_use = "theme_blank")
macrophage
ggsave("3.macrophagedimplot.tiff",macrophage,width = 6, height = 4, dpi =500)


endolithial <- FeatureDimPlot(Srt.GSE227075, features = c("Pecam1","Eng","Cdh5"), theme_use = "theme_blank")
endolithial
ggsave("4.endolithialdimplot.tiff",endolithial,width = 6, height = 4,dpi =500)


Celltypeheatmap <- GroupHeatmap(
  srt = Srt.GSE227075,
  features = c(
    "Krt10","Krt1", #Epi
    "Cd3d","Cd3e","Cd3g",
    "Hbb-bt", "Hba-a1","Hba-a2",#Redbloodcell
    "Dcn",
    "Postn",
    "Myog",
    "Des",
    "Myod1",
    "Pecam1","Eng","Cdh5",#endolithial
    "Arg1","Msr1","Mrc1",#Mac
    "Thbs4",
    "Fmod",
    "Lox",
    "Ctsk",
    "Acp5",
    "S100a9",
    "S100a8",
    "Cd74",
    "Cd79a",
    "Hmgb2",
    "Stmn1",
    "Top2a",
    "Gldn",
    "Cryab",
    "Rgs5",
    "Myl9",
    "Ccr9",
    "Siglech",
    "Klk1",
    "Adipoq",
    "Plin1",
    "Car3"
  ),
  group.by = c("Celltype"),
  flip = TRUE,
  cell_annotation = c("Batchid"),
  heatmap_palette = "YlOrRd",
  show_row_names = FALSE, row_names_side = "left",
  show_column_names = TRUE, column_names_side = "bottom",
  column_names_rot = 45,
  add_dot = TRUE, add_bg = TRUE,
  dot_size = unit(5, "mm"),
  width = 8,height = 6,
  nlabel = 0,
)
Celltypeheatmap
ggsave("Celltypeheatmap.pdf",Celltypeheatmap$plot,width = 15, height = 8)


#macrophage subset
Srt.GSE227075 <- SetIdent(Srt.GSE227075, value = 'Celltype')
Srt.macrophage <- subset(Srt.GSE227075, Celltype == "Macrophages")
table(Srt.macrophage@meta.data$Celltype)
DimPlot(object = Srt.macrophage, reduction = 'umap', group.by = "Batchid", pt.size=0.1)
Srt.macrophage.recluster <- RunPCA(object = Srt.macrophage, verbose = FALSE)
ElbowPlot(object = Srt.macrophage.recluster, ndims = 50)
use.pcs = 1:20
Srt.macrophage.recluster <- FindNeighbors(object = Srt.macrophage.recluster, reduction = "harmony", dims = use.pcs)
Srt.macrophage.recluster <- FindClusters(object = Srt.macrophage.recluster, reduction.type = "harmony", dims.use = use.pcs, resolution = seq(0.1,2,0.1), print.output = FALSE, save.SNN = TRUE)
Srt.macrophage.recluster <- RunUMAP(object = Srt.macrophage.recluster, reduction = "harmony", dims = use.pcs)
sapply(grep("^RNA_snn_res",colnames(Srt.macrophage.recluster@meta.data),value = TRUE), function(x) length(unique(Srt.macrophage.recluster@meta.data[,x])))
p3 <- CellDimPlot(srt = Srt.macrophage.recluster, group.by = "RNA_snn_res.2", reduction = "UMAP",label = TRUE,theme_use = "theme_blank")
p3
FeatureDimPlot(Srt.macrophage.recluster, features = c("Hdc"), theme_use = "theme_blank")


CellDimPlot(srt = Srt.macrophage.recluster, group.by = "Batchid", reduction = "UMAP", theme_use = "theme_blank")
ggsave("macrophagecluster.pdf",p3,width = 6, height = 4)
DimPlot(object = Srt.macrophage.recluster, reduction = 'umap', group.by = "RNA_snn_res.2", pt.size=0.1, label = TRUE, repel = TRUE)
DimPlot(object = Srt.macrophage.recluster, reduction = 'umap', group.by = "RNA_snn_res.2", split.by = "Batchid",  pt.size=0.1, label = TRUE, ncol = 2)
Srt.macrophage.recluster <- SetIdent(Srt.macrophage.recluster, value = 'RNA_snn_res.2')
plan("multisession", workers = 6)
markers_all_RNA <- FindAllMarkers(object = Srt.macrophage.recluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "LR")
write.csv(markers_all_RNA, file = "macrophage_subset_res2.csv")
plan("sequential")


#cell ratio barplot
Cellratio.macrophage <- prop.table(table(Idents(Srt.macrophage.recluster), Srt.macrophage.recluster$Batchid), margin = 2)
Cellratio.macrophage <- as.data.frame(Cellratio.macrophage)
colourCount = length(unique(Cellratio.macrophage$Var1))
levels(Cellratio.macrophage$Var1)
Cellratio.macrophage$Var2 <- factor(Cellratio.macrophage$Var2,levels = c("Sham", "Day1", "Day3", "Day7"))
library(ggplot2)
partition.barplot <- ggplot(Cellratio.macrophage, aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.8,size = 0.5,colour = '#222222') + 
  geom_col(size = 0.5)+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_y_reverse()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  geom_text(aes(label = Var1), 
            position = position_stack(vjust = 0.5))
partition.barplot


#Subtype annotation
table(Srt.macrophage.recluster$RNA_snn_res.2)
Srt.macrophage.recluster <- SetIdent(Srt.macrophage.recluster, value = 'RNA_snn_res.2')
Srt.macrophage.recluster$Subtype <- NULL
new.cluster.ids <- c("NULL",#0
                     "Trem2+",
                     "Trem2+",
                     "Trem2+",
                     "Mrc1+",
                     "Mrc1+",
                     "Mrc1+",
                     "Trem2+",
                     "Cd52+",
                     "Cd52+",
                     "Arg1+",#10
                     "Cd52+",
                     "Trem2+",
                     "Mrc1+",
                     "Mrc1+",
                     "NULL",#15
                     "Cd52+",
                     "Hdc+",
                     "Cd52+",
                     "Trem2+",
                     "Cd52+",
                     "Arg1+",
                     "NULL",#22
                     "Trem2+",
                     "Trem2+",
                     "Mrc1+",
                     "Trem2+",
                     "Mrc1+",
                     "NULL",#28
                     "Mrc1+")
names(new.cluster.ids) <- levels(Srt.macrophage.recluster) 
Srt.macrophage.recluster <- RenameIdents(Srt.macrophage.recluster, new.cluster.ids)
Srt.macrophage.recluster[["Subtype"]] <- Idents(Srt.macrophage.recluster)
table(Srt.macrophage.recluster$Subtype)

#filter null
Srt.macrophage.recluster <- subset(Srt.macrophage.recluster, Subtype %in% c("Trem2+", "Mrc1+", "Cd52+", "Arg1+", "Hdc+"))
table(Srt.macrophage.recluster$Subtype)
annotation <- list(
  "Trem2+" = c('NULL','Trem2+'),
  'Mrc1+' = 'Mrc1+',
  'Cd52+' = 'Cd52+',
  'Arg1+' = 'Arg1+',
  'Hdc+' = 'Hdc+')
Srt.macrophage.recluster <- RenameClusters(Srt.macrophage.recluster, group.by = "Subtype", nameslist = annotation, name = "Subtype")
table(Srt.macrophage.recluster$Subtype)
p4 <- CellDimPlot(srt = Srt.macrophage.recluster, group.by = "Subtype", reduction = "UMAP", theme_use = "theme_blank")
p4
ggsave("macrophageclusterSubtype.pdf",p4,width = 6, height = 4)


Srt.macrophage.recluster <- SetIdent(Srt.macrophage.recluster, value = 'Subtype')
p5 <- GroupHeatmap(
  srt = Srt.macrophage.recluster,
  features = c(
    "Trem2",
    "Mrc1",
    "Cd52",
    "Arg1",
    "Hdc", "Vegfa", "Mmp9", "Cxcr2","Il1b"
  ),
  group.by = c("Subtype"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase", "G2M.Score"),
  show_row_names = FALSE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE
)
print(p5$plot)
ggsave("Groupheatmap.pdf",p5$plot,width = 6, height = 6)



features <- c("Hdc")
p6 <- FeatureDimPlot(Srt.macrophage.recluster, features = features, theme_use = "theme_blank")
p6
ggsave("HdcCxcr2inmacrophage.pdf",p6,width = 6, height = 4)


features <- c("Hdc","Cxcr2","Vegfa","Mmp9","Il1b")
p7 <- FeatureDimPlot(Srt.macrophage.recluster, features = features, theme_use = "theme_blank")
p7
ggsave("Hdcgenesinmacrophage.pdf",p7,width = 6, height = 4)


features <- c("Hdc","Cxcr2","Vegfa","Mmp9")
FeaturePlot(Srt.macrophage.recluster, features = features)


p8 <- FeatureDimPlot(
  srt = Srt.macrophage.recluster, features = c("Cxcr2", "Hdc","Vegfa","Mmp9","Il1b"),
  compare_features = TRUE, label = FALSE, label_insitu = TRUE,
  reduction = "UMAP", theme_use = "theme_blank"
)
p8
ggsave("Hdcgenesinmacrophagetotal.pdf",p8,width = 6, height = 4)

#DEGtest and VolcanoPlot
register(SnowParam(workers = 2, progressbar = TRUE))
Srt.macrophage.recluster <- RunDEtest(srt = Srt.macrophage.recluster, group_by = "Subtype", fc.threshold = 1, only.pos = FALSE)
p9 <- VolcanoPlot(srt = Srt.macrophage.recluster, group_by = "Subtype", features_label = c("Hdc"))
p9
ggsave("Volcanoplot.pdf",p9,dpi = 1000)


#big heatmap
DEGs <- Srt.macrophage.recluster@tools$DEtest_Subtype$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
# Annotate features with transcription factors and surface proteins
Srt.macrophage.recluster <- AnnotateFeatures(Srt.macrophage.recluster, species = "Mus_musculus", db = c("TF", "CSPA"))
p10 <- FeatureHeatmap(
  srt = Srt.macrophage.recluster, group.by = "Subtype", features = DEGs$gene, feature_split = DEGs$group1,
  max_cells = 200,
  species = "Mus_musculus", db = c("GO_BP", "KEGG", "WikiPathway"), anno_terms = TRUE,
  feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")))
print(p10$plot)
ggsave("bigheatmap.pdf",p10$plot,limitsize = FALSE, width = 20, height = 8)

#cellratio barplot
p11 <- CellStatPlot(srt = Srt.macrophage.recluster, stat.by = "Subtype", group.by = "Batchid", label = FALSE,
                    stat_type = "percent")
p11
ggsave("cellratiobarplot.pdf",p11,dpi = 1000)

p12 <- CellStatPlot(Srt.macrophage.recluster, stat.by = c("Subtype", "Batchid"), plot_type = "chord")
p12
ggsave("cellcountchordplot.pdf",p12,dpi = 1000)


#pathway
library(ggkegg)
g <- pathway("mmu04060") |> mutate(marker_1=append_cp(ekegg))
gg <- ggraph(g, layout="manual", x=x, y=y)+
  geom_node_rect(aes(filter=marker_1), fill="tomato")+ ## Marker 1
  overlay_raw_map("mmu04060", transparent_colors = c("#cccccc","#FFFFFF","#BFBFFF","#BFFFBF"))+
  theme_void()
gg
ggsave(filename = "cytokine.pdf", plot = gg, width = 9, height = 9, dpi = 300)


#CytoTRACE
plan("sequential")
cytotrace.matrix <- as.matrix(Srt.macrophage.recluster@assays$RNA@counts)
cytotrace.matrix <- cytotrace.matrix[apply(cytotrace.matrix > 0,1,sum) >= 5,]
cytotrace.results <- CytoTRACE(cytotrace.matrix,ncores = 1, enableFast = FALSE)
cytotrace.clusters <- Srt.macrophage.recluster$Subtype
cytotrace.clusters <- as.character(cytotrace.clusters)
names(cytotrace.clusters) <- rownames(Srt.macrophage.recluster@meta.data)
emb <- Srt.macrophage.recluster@reductions[["umap"]]@cell.embeddings
plotCytoTRACE(cytotrace.results, phenotype = cytotrace.clusters, emb = emb, outputDir = './')
plotCytoGenes(cytotrace.results, numOfGenes = 30, outputDir = './')


#monocle
CellDimPlot(Srt.macrophage.recluster, group.by = "Batchorder", reduction = "UMAP", label = TRUE, theme_use = "theme_blank")
Srt.macrophage.recluster <- RunMonocle3(srt = Srt.macrophage.recluster, reduction = "umap")
names(Srt.macrophage.recluster@tools$Monocle3)
trajectory <- Srt.macrophage.recluster@tools$Monocle3$trajectory
milestones <- Srt.macrophage.recluster@tools$Monocle3$milestones

CellDimPlot(Srt.macrophage.recluster, group.by = "Subtype", reduction = "UMAP", label = TRUE, theme_use = "theme_blank") + trajectory + milestones
CellDimPlot(Srt.macrophage.recluster, group.by = "Monocle3_clusters", reduction = "UMAP", label = TRUE, theme_use = "theme_blank") + trajectory
FeatureDimPlot(Srt.macrophage.recluster, features = "Monocle3_Pseudotime", reduction = "UMAP", theme_use = "theme_blank") + trajectory


CellDimPlot(Srt.macrophage.recluster, group.by = "partitions", reduction = "UMAP", label = TRUE, theme_use = "theme_blank") + trajectory + milestones

#L1
cds <- Srt.macrophage.recluster@tools$Monocle3$cds
cds_sub1 <- monocle3::choose_graph_segments(cds, starting_pr_node = 34, ending_pr_nodes = 485)
Srt.macrophage.recluster$Lineages_1 <- NA
Srt.macrophage.recluster$Lineages_1[colnames(cds_sub1)] <- Srt.macrophage.recluster$Monocle3_Pseudotime[colnames(cds_sub1)]
#L2
cds <- Srt.macrophage.recluster@tools$Monocle3$cds
cds_sub2 <- monocle3::choose_graph_segments(cds, starting_pr_node = 34, ending_pr_nodes = 255)
Srt.macrophage.recluster$Lineages_2 <- NA
Srt.macrophage.recluster$Lineages_2[colnames(cds_sub2)] <- Srt.macrophage.recluster$Monocle3_Pseudotime[colnames(cds_sub2)]
#L3
cds <- Srt.macrophage.recluster@tools$Monocle3$cds
cds_sub3 <- monocle3::choose_graph_segments(cds, starting_pr_node = 34, ending_pr_nodes = 83)
Srt.macrophage.recluster$Lineages_3 <- NA
Srt.macrophage.recluster$Lineages_3[colnames(cds_sub3)] <- Srt.macrophage.recluster$Monocle3_Pseudotime[colnames(cds_sub3)]
table(Srt.macrophage.recluster$Lineages_3)

p13 <- CellDimPlot(Srt.macrophage.recluster, group.by = "Subtype", lineages = c("Lineages_1","Lineages_2","Lineages_3"), theme_use = "theme_blank")
p13
ggsave(filename = "pseudotime.pdf", plot = p13, width = 6, height = 6, dpi = 300)


#WCGNA
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
Track_genes <- Track_genes[,c(2,3,4,5)] %>% filter(Track_genes$q_value < 1e-3)
write.csv(Track_genes, "Trajectory_genes.csv", row.names = F)
write.csv(Track_genes, "Trajectory_genes_name.csv", row.names = T)


p14 <- plot_genes_in_pseudotime(cds[features,], color_cells_by="Subtype",min_expr=0.5, ncol = 2)
p14
ggsave(filename = "hdcgenesinpseudotime.pdf", plot = p14, width = 8, height = 6, dpi = 1000)

gene_short_name = rownames(cds)
cds <- preprocess_cds(cds)
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 6)
write.csv(gene_module, "Genes_Module.csv", row.names = F)
cell_group <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$Subtype)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
p15 <- pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
p15
ggsave("Genes_Module.pdf", plot = p15, width = 10, height = 10)

#save
saveRDS(Srt.macrophage.recluster,"smr.rds")
Srt.macrophage.recluster <- readRDS("smr.rds")
Srt.GSE227075 <- readRDS("srtGSE227075.rds")

#pool subtype into the total
subtype_levels <- levels(Srt.macrophage.recluster$Subtype)
celltype_levels <- levels(Srt.GSE227075$Celltype)
new_levels <- unique(c(celltype_levels, subtype_levels))
levels(Srt.GSE227075$Celltype) <- new_levels
Srt.GSE227075$Celltype[match(colnames(Srt.macrophage.recluster), colnames(Srt.GSE227075))] <- Srt.macrophage.recluster$Subtype
table(Srt.GSE227075$Celltype)


Srt.GSE227075 <- SetIdent(Srt.GSE227075, value = 'Celltype')
Srt.GSE227075.cellchat <- subset(Srt.GSE227075, idents = c("FAPs",
                                                  "B cells",
                                                  "MuSCs",
                                                  "Neutrophils",
                                                  "Cycling basal cells",
                                                  "Enodothelial cells",
                                                  "Schwann cells",
                                                  "NK/T cells",
                                                  "Tenocytes",
                                                  "Erythrocytes",
                                                  "Epithelial cells",
                                                  "Trem2+",
                                                  "Mrc1+",
                                                  "Cd52+",
                                                  "Arg1+",
                                                  "Hdc+"))
table(Srt.GSE227075.cellchat$Celltype)

new_idents <- Idents(Srt.GSE227075)
new_idents[new_idents == "Macrophages"] <- NULL
Idents(Srt.GSE227075) <- new_idents

table(Srt.GSE227075$Celltype)



#cellchat
Srt.GSE227075.cellchat <- NormalizeData(Srt.GSE227075.cellchat, verbose = F)
Srt.GSE227075.cellchat <- ScaleData(Srt.GSE227075.cellchat)
Srt.GSE227075.cellchat <- SetIdent(Srt.GSE227075.cellchat, value = 'Celltype')

data.input <- GetAssayData(Srt.GSE227075.cellchat, assay = "RNA", slot = "data")
labels <- Idents(Srt.GSE227075.cellchat)
meta <- data.frame(group = labels, row.names = names(labels))
cellChat <- createCellChat(object = Srt.GSE227075.cellchat, group.by = "ident", assay = "RNA")

CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellChat@DB <- CellChatDB.use

cellChat <- subsetData(cellChat)
future::plan("multisession", workers = 4)
cellChat <- identifyOverExpressedGenes(cellChat)
cellChat <- identifyOverExpressedInteractions(cellChat)
cellChat <- computeCommunProb(cellChat, type = "triMean")
cellChat <- filterCommunication(cellChat, min.cells = 10)
cellChat <- computeCommunProbPathway(cellChat)
cellChat <- aggregateNet(cellChat)

pathways.show.all <- cellChat@netP$pathways
table(pathways.show.all)
# select one pathway
pathways.show <- c("VEGF") 
par(mfrow=c(1,1))

pdf("vegfnetchord.pdf")
vegfnet <- netVisual_aggregate(cellChat, signaling = pathways.show, layout = "chord",color.use = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL)
dev.off()

ggsave(filename = "vegfnet.tiff", plot = vegfnet, dpi = 300)



table(cellChat@meta$labels)
netVisual_bubble(cellChat, sources.use = 17, targets.use = c(1,4,6:12), remove.isolate = FALSE)


#NIchenet
library(circlize)
library(nichenetr)


ligand_target_matrix = readRDS("C:/Users/zjf/Downloads/ligand_target_matrix_nsga2r_final_mouse.rds")
lr_network = readRDS("C:/Users/zjf/Downloads/lr_network_mouse_21122021.rds")
weighted_networks = readRDS("C:/Users/zjf/Downloads/weighted_networks_nsga2r_final_mouse.rds")
lr_network = lr_network %>% distinct(from, to)
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))


table(Srt.GSE227075$Celltype)
Srt.GSE227075 <- SetIdent(Srt.GSE227075, value = 'Celltype')

receiver = "Enodothelial cells"
expressed_genes_receiver = get_expressed_genes(receiver, Srt.GSE227075, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

sender_celltypes = c("Trem2+","Mrc1+", "Cd52+", "Arg1+", "Hdc+")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, Srt.GSE227075, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

seurat_obj_receiver= subset(Srt.GSE227075, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["Batchid", drop=TRUE]])
condition_oi = c("Day1","Day3")
condition_reference = "Sham" 

DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]



ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()


ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

best_upstream_ligands = ligand_activities %>% top_n(60, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
nichenet <- DotPlot(Srt.macrophage.recluster, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
ggsave("nichenet.tiff",nichenet, dpi = 300, width = 16, height = 6)




specified_ligands <- c("Vegfa", "Mmp9", "Il1b", "Il1a", "Il15", "Il16")
best_upstream_ligands = ligand_activities %>%
  filter(test_ligand %in% specified_ligands | rank(-aupr_corrected) <=20) %>%
  arrange(-aupr_corrected) %>%
  pull(test_ligand) %>%
  unique()


avg_expression_ligands = AverageExpression(Srt.macrophage.recluster, features = best_upstream_ligands)
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()

sender_ligand_assignment[1:5,1:5]

sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)

(sender_ligand_assignment)

all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = best_upstream_ligands %>% setdiff(unique_ligands)



Trem2_specific_ligands = sender_ligand_assignment$`Trem2+` %>% names() %>% setdiff(general_ligands)
Cd52_specific_ligands = sender_ligand_assignment$`Cd52+` %>% names() %>% setdiff(general_ligands)
Arg1_specific_ligands = sender_ligand_assignment$`Arg1+` %>% names() %>% setdiff(general_ligands)
Hdc_specific_ligands = sender_ligand_assignment$`Hdc+` %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("Trem2-specific", times = Trem2_specific_ligands %>% length()),
                  rep("Cd52-specific", times = Cd52_specific_ligands %>% length()),
                  rep("Arg1-specific", times = Arg1_specific_ligands %>% length()),
                  rep("Hdc-specific", times = Hdc_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  ligand = c(Trem2_specific_ligands, Cd52_specific_ligands, Arg1_specific_ligands, Hdc_specific_ligands, general_ligands))

ligand_type_indication_df %>% head



active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "geneset_oi") %>% inner_join(ligand_type_indication_df)
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.78)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

circos_links



grid_col_ligand =c("General" = "green4",
                   "Trem2-specific" = "lightgreen",
                   "Cd52-specific" = "violet",
                   "Arg1-specific" = "steelblue2",
                   "Hdc-specific" = "royalblue")
grid_col_target =c(
  "geneset_oi" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% dplyr::select(ligand, target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)



target_order = circos_links$target %>% unique()
ligand_order = c(Trem2_specific_ligands, Cd52_specific_ligands, Arg1_specific_ligands, Hdc_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)

transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency


width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Trem2-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Cd52-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Arg1-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Hdc-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_ligand_target,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "geneset_oi") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)

tiff("ligand_target_circos.tiff", width = 10, height = 10, units = 'in', res = 300)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.3)
}, bg.border = NA) #
circos.clear()
dev.off()



#Endolithial subset
Srt.GSE227075 <- SetIdent(Srt.GSE227075, value = 'Celltype')
table(Srt.GSE227075$Celltype)
Srt.endolithial <- subset(Srt.GSE227075, Celltype == "Enodothelial cells")

Srt.endolithial <- SetIdent(Srt.endolithial, value = 'Batchid')
original_idents <- Idents(object = Srt.endolithial)
Condition <- ifelse(original_idents %in% c("Day1", "Day3", "Day7"), "LID", 
                     ifelse(original_idents == "Sham", "Sham", NA))
Condition[is.na(Condition)] <- "Other"
Idents(object = Srt.endolithial) <- Condition
table(Idents(Srt.endolithial))

Srt.endolithial[["Condition"]] <- Idents(Srt.endolithial)




Srt.endolithial <- RunDEtest(srt = Srt.endolithial, group_by = "Condition", fc.threshold = 1, only.pos = FALSE)
p15 <- VolcanoPlot(srt = Srt.endolithial, group_by = "Condition", features_label = c("Cxcl2"),x_metric = "avg_log2FC")
p15
ggsave("Cxcl2Volcanoplot.pdf",p15,dpi = 1000)

DEGs.endothelial <- Srt.endolithial@tools$DEtest_Condition$AllMarkers_wilcox
write.csv(DEGs.endothelial,"DEGs.endothelial.csv")


