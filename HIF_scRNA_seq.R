# Load the various packages required
library(ReactomePA)
library(ComplexHeatmap)
library(clusterProfiler)
library(circlize)
library(org.Mm.eg.db)
library(tidyverse)
library(Seurat)
library(flextable)
library(officer)
library(plyr)

#Creation of the Seurat object WT 5 months
data2_10x = Read10X(data.dir = " ")
Ech1 = CreateSeuratObject(data2_10x,
                          project = "WT",
                          assay = "RNA",
                          min.cells = 10,
                          min.features = 100,
                          names.field = 1,
                          names.delim = "_",
                          meta.data = NULL
)
Ech1$orig.ident="1_WT"
Ech1[["percent.mt"]] <- PercentageFeatureSet(Ech1, pattern = "^mt-")
VlnPlot(Ech1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=1, ncol = 3)
ggsave("VlnPlot_WT.tiff")

#Creation of the Seurat object PTEN_3mo after gene invalidation
data2_10x = Read10X(data.dir = " ")
Ech2 = CreateSeuratObject(data2_10x,
  project = "3moPTEN",
  assay = "RNA",
  min.cells = 10,
  min.features = 100,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL
)
Ech2$orig.ident="3mo"
Ech2[["percent.mt"]] <- PercentageFeatureSet(Ech2, pattern = "^mt-")
VlnPlot(Ech2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=1, ncol = 3)
ggsave("VlnPlot_3mo.tiff")
remove(data2_10x)

#Creation of the combined Seurat object 
Ech.list <- list(Ech1, Ech2)
PTEN_3mo <- lapply(X = Ech.list, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 50)
  x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
remove(Ech1, Ech2, Ech.list)

PTEN_3mo <- FindIntegrationAnchors(object.list = PTEN_3mo, dims = 1:30)
PTEN_3mo <- IntegrateData(anchorset = PTEN_3mo, dims = 1:30)


#Visualisation du tSNE
DefaultAssay(PTEN_3mo) <- "integrated"

# Run the standard workflow for visualization and clustering
PTEN_3mo <- ScaleData(PTEN_3mo, verbose = FALSE)
PTEN_3mo <- RunPCA(PTEN_3mo, npcs = 50, verbose = T, ndims.print = 1:15 )

# t-SNE and Clustering

PTEN_3mo <- FindNeighbors(PTEN_3mo, reduction = "pca", dims = 1:30)
PTEN_3mo <- FindClusters(PTEN_3mo,resolution =  0.6)
PTEN_3mo <- RunTSNE(PTEN_3mo, dims.use = 1:30, reduction.use = "pca", perplexity = 30)
PTEN_3mo <- RunUMAP(PTEN_3mo, dims = 1:30)
PTEN_3mo$CellType <- Idents(PTEN_3mo)

#distribution of the various clusters
table <-(data.frame(table(PTEN_3mo@meta.data$CellType, PTEN_3mo@meta.data$orig.ident)))
table <- table %>% pivot_wider(names_from = "Var2", values_from="Freq")
colnames(table)[1]<-"cluster"
table$'1_WT_percent' <- round(100 * (table$'1_WT'/sum(table$'1_WT')),1)
table$'3mo_percent' <- round(100 * (table$'3mo'/sum(table$'3mo')),1)
x <- sum(table$'3mo')+sum(table$'1_WT')
table$'percent_total' <- round(100 * ((table$'3mo'+table$'1_WT')/x),1)
write.table(table, "tavle_cluster_WTPTEN.txt", row.names = T)
table_PTEN_WT <- table
remove(x, table)

#elimination of the clusters containing less than 0.5% of the cells (i.e. 25, 26)
PTEN_3mo <- subset(PTEN_3mo, idents=c(25, 26), invert=T)

#annotation of the clusters
#Supp Figure 1A
DimPlot(PTEN_3mo, reduction = "umap", label = T,  label.size = 3, pt.size = 0.05, split.by="orig.ident") + theme(axis.text = element_text( size = 8), axis.title = element_text(size = 8), legend.text=element_text( size = 8),title = element_text( size = 8)) + NoLegend()
ggsave("Supp_Fig1A.eps", units = "mm", height = 80, width = 80, dpi=300)

#Supp Figure 1B
SuppFigure1B <- c("Epcam", "Ptprc", "Vim")
DotPlot(dot.min = 0.2, PTEN_3mo, features = SuppFigure1B, assay = "RNA") + theme(axis.text.x = element_text(size = 7, hjust = 1, vjust = 1, angle=45), axis.text.y = element_text(size = 7), axis.title = element_text(size = 8), title = element_text( size = 7), axis.title.x = element_blank(), legend.title = element_text(size = 7), legend.key.height = unit(0.4, "lines"), legend.key.width = unit(0.5, "lines"), legend.text = element_text(size = 7)) + guides(size = guide_legend(title = 'Percent'), color = guide_colorbar(title = 'Expression')) 
ggsave("SuppFig1B.eps", units = "mm", height = 90, width = 50, dpi=300)

#novel annotation 
new.ids <- c("stroma", "stroma", "stroma", "stroma", "leukocytes", 
             "Luminal-A", "stroma", "stroma", "leukocytes", "T-cells", 
             "leukocytes", "Luminal-C", "Luminal-C", "endo", "Basal", 
             "Luminal-B",  "stroma", "leukocytes", "Luminal-A", "endo",
             "Luminal-A", "leukocytes", "MDSC", "B-cells", "Luminal-B")
names(new.ids) <- levels(PTEN_3mo)
PTEN_3mo <- RenameIdents(PTEN_3mo, new.ids)
my_levels <- rev(c("Luminal-A",  "Luminal-B",  "Luminal-C", "Basal", "leukocytes" , "stroma","endo", "MDSC",     "B-cells", "T-cells"))
PTEN_3mo@active.ident <- factor(PTEN_3mo@active.ident, levels=my_levels)
PTEN_3mo$fig1 <- PTEN_3mo@active.ident
Idents(PTEN_3mo) <- "fig1"

#Figure 1A
DimPlot(PTEN_3mo, reduction = "umap", label = F, pt.size = 0.1, split.by="orig.ident", order = rev(c("Luminal-A",  "Luminal-B",  "Luminal-C","Basal", "leukocytes" , "stroma","endo", "MDSC",     "B-cells", "T-cells"))) + theme(axis.text = element_text( size = 8), axis.title = element_text(size = 8), legend.text=element_text( size = 8),title = element_text( size = 8))
ggsave("Fig1A.eps", units = "mm", height = 60, width = 100, dpi=300)

#Figure 1B
Figure1B <- c("Krt5","Krt8", "Nkx3-1","Tmprss2", "Pbsn", "Ly6a", "Tacstd2","Pate4", "Krt4", "Pecam1", "S100a8", "Cd3e", "Cd79a", "Col1a1", "Adgre1", "Itgam")
DotPlot(dot.min = 0.2, PTEN_3mo, features = Figure1B, assay = "RNA") + theme(axis.text.x = element_text(size = 7, hjust = 1, vjust = 0.5, angle=90), axis.text.y = element_text(size = 7), axis.title = element_text(size = 8), title = element_text( size = 7), axis.title.x = element_blank(), legend.title = element_text(size = 7), legend.key.height = unit(0.4, "lines"), legend.key.width = unit(0.5, "lines"), legend.text = element_text(size = 7)) + guides(size = guide_legend(title = 'Percent'), color = guide_colorbar(title = 'Expression')) 
ggsave("Fig1B.eps", units = "mm", height = 60, width = 80, dpi=300)

#Table_markers
Idents(PTEN_3mo) <- 'CellType'
Signatures_PTEN<-data.frame()

for (i in 1:length(levels(PTEN_3mo))) 
  {
z <- FindConservedMarkers(PTEN_3mo,ident.1 = sort(levels(PTEN_3mo))[i], grouping.var="orig.ident", logfc.threshold = 0.25, min.pct = 0.5, verbose = T, only.pos = T, assay = 'integrated')
z$gene <-rownames(z)
z$NCBI_id <- mapIds(org.Mm.eg.db, z$gene, 'ENTREZID', 'SYMBOL')
z$cluster <- sort(levels(PTEN_3mo))[[i]]
Signatures_PTEN<- rbind(Signatures_PTEN, z)
}
Supplementary_Table_1 <- Signatures_PTEN
write.csv(x=Signatures_PTEN, file = "Supplementary_Table_1.csv", row.names = F) 

#Determine DEG PTEN vs. WT for luminal C and A
Idents(PTEN_3mo) <- 'fig1'
z = FindMarkers(PTEN_3mo, ident.1 = '3mo', ident.2 = '1_WT', group.by='orig.ident', subset.ident=c('Luminal-A', 'Luminal-C'), assay='RNA', verbose = TRUE)
  z$gene <-rownames(z)
  z$NCBI_id <- mapIds(org.Mm.eg.db, z$gene, 'ENTREZID', 'SYMBOL')
  z <- z %>% relocate(gene)
  write.csv(x=z, file = "Supplementary_Table_2.csv", row.names = F) 
Supplementary_Table_2 = z
  Gene_up <- z[z$p_val<0.05 & z$avg_log2FC>0,]   
  KEGG_luminal <- enrichKEGG(Gene_up$NCBI_id, organism = "mmu", keyType="kegg", pvalueCutoff=0.05, pAdjustMethod = "fdr", universe,  minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.05)
KEGG_luminal<-setReadable(KEGG_luminal, org.Mm.eg.db, keyType = "ENTREZID")

my_gene_list <- KEGG_luminal@result[KEGG_luminal@result$Description=='HIF-1 signaling pathway',]
my_gene_list <- str_split(my_gene_list$geneID, "/", simplify=T)
Supplementary_Table_3 <- KEGG_luminal@result
write.csv(x = KEGG_luminal@result, file = "Supplementary_Table_3.csv", sep = " ", row.names = F)

#Figure 1D 
Idents(PTEN_3mo) <- "fig1"
new.ids <- c("T-cells",    "B-cells",    "MDSC",       "endo",       "stroma",     "leukocytes", "Basal"    , "Luminal-A/C",  "Luminal-B",  "Luminal-A/C")
names(new.ids) <- levels(PTEN_3mo)
PTEN_3mo <- RenameIdents(PTEN_3mo, new.ids)
PTEN_3mo$fig2 <- PTEN_3mo@active.ident
Idents(PTEN_3mo) <- "fig2"

for (i in 1:(length(my_gene_list))) {
  
       VlnPlot(PTEN_3mo, features = my_gene_list[[i]],pt.size = 0.1, assay = "RNA", idents="Luminal-A/C", split.by = "orig.ident") + theme(axis.title = element_text(size = 8), axis.line.x = element_line(size = 0), axis.line.y.left = element_line(size = 0.2), axis.text.y.left =  element_text(size = 8), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_line(size=0.2), axis.title.y = element_blank(), axis.text.x = element_text(size = 8, angle=0, hjust=0.5, vjust=0), title = element_text(size = 8)) + NoLegend()
name <- paste0(my_gene_list[[i]],"_Fig1D.eps")
ggsave(name, units = "mm", height = 40, width = 40, dpi=300)  
}

#Creation of the Seurat object PTEN/HIF_3mo after gene invalidation
data2_10x = Read10X(data.dir = "~/Documents/SingleCell/outs/MNGC246_filtered_feature_bc_matrix/")
Ech3 = CreateSeuratObject(data2_10x,
  project = "3moHIF",
  assay = "RNA",
  min.cells = 10,
  min.features = 100,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL
)
Ech3$orig.ident="3moHIF"
Ech3[["percent.mt"]] <- PercentageFeatureSet(Ech3, pattern = "^mt-")
VlnPlot(Ech3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=1, ncol = 3)
ggsave("VlnPlot_HIF3mo.tiff")
remove(data2_10x)

#Creation of the Seurat object PTEN_3mo after gene invalidation
data2_10x = Read10X(data.dir = "~/Documents/SingleCell/outs/MNGC244_filtered_feature_bc_matrix/")
Ech2 = CreateSeuratObject(data2_10x,
  project = "3moPTEN",
  assay = "RNA",
  min.cells = 10,
  min.features = 100,
  names.field = 1,
  names.delim = "_",
  meta.data = NULL
)
Ech2$orig.ident="3mo"
Ech2[["percent.mt"]] <- PercentageFeatureSet(Ech2, pattern = "^mt-")
VlnPlot(Ech2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=1, ncol = 3)
ggsave("VlnPlot_3mo.tiff")
remove(data2_10x)

#Creation of the combined Seurat object 
Ech.list <- list(Ech2, Ech3)
HIF <- lapply(X = Ech.list, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 50)
  x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
remove(Ech3, Ech2, Ech.list)

HIF <- FindIntegrationAnchors(object.list = HIF, dims = 1:30)
HIF <- IntegrateData(anchorset = HIF, dims = 1:30)


#Visualisation du tSNE
DefaultAssay(HIF) <- "integrated"

# Run the standard workflow for visualization and clustering
HIF <- ScaleData(HIF, verbose = FALSE)
HIF <- RunPCA(HIF, npcs = 50, verbose = T, ndims.print = 1:15 )

# t-SNE and Clustering

HIF <- FindNeighbors(HIF, reduction = "pca", dims = 1:30)
HIF <- FindClusters(HIF,resolution =  0.6)
HIF <- RunTSNE(HIF, dims.use = 1:30, reduction.use = "pca", perplexity = 30)
HIF <- RunUMAP(HIF, dims = 1:30)

#distribution of the various clusters
table <-(data.frame(table(HIF@meta.data$CellType, HIF@meta.data$orig.ident)))
table <- table %>% pivot_wider(names_from = "Var2", values_from="Freq")
colnames(table)[1]<-"cluster"
table$'1_WT_percent' <- round(100 * (table$'3mo'/sum(table$'3mo')),1)
table$'3mo_percent' <- round(100 * (table$'3moHIF'/sum(table$'3moHIF')),1)
x <- sum(table$'3mo')+sum(table$'3moHIF')
table$'percent_total' <- round(100 * ((table$'3mo'+table$'3moHIF')/x),1)
table_HIF <- table
write.table(table, "Table_HIF.txt", row.names = T)
remove(x, table)

#elimination of the clusters containing less than 0.5% of the cells (i.e. 26, 27, 27, 29)
HIF <- subset(HIF, idents=c(26, 27, 28, 29), invert=T)

#novel annotation 
HIF$CellType <- HIF@active.ident
new.ids <- c("stroma", "leukocytes", "stroma", "stroma", "Basal", 
             "stroma", "endo", "Luminal-A", "MDSC", "leukocytes", 
             "stroma", "leukocytes", "Luminal-B", "T-cells", "Luminal-C", 
             "leukocytes",  "Luminal-C", "T-cells", "T-cells", "T-cells",
             "B-cells", "leukocytes", "Luminal-A", "endo", "stroma", "Luminal-C")
names(new.ids) <- levels(HIF)
HIF <- RenameIdents(HIF, new.ids)
my_levels <- c("Luminal-A",  "Luminal-B",  "Luminal-C","Basal", "leukocytes" , "stroma","endo", "MDSC",     "B-cells", "T-cells")
HIF@active.ident <- factor(HIF@active.ident, levels=my_levels)
HIF$fig1 <- HIF@active.ident
Idents(HIF) <- "fig1"

#Figure 1A
Idents(HIF) <- 'fig1'
DimPlot(HIF, reduction = "umap", label = F,  label.size = 2, pt.size = 0.1, split.by="orig.ident") + theme(axis.text = element_text( size = 8), axis.title = element_text(size = 8), legend.text=element_text( size = 8),title = element_text( size = 8))
ggsave("Fig3A.eps", units = "mm", height = 60, width = 100, dpi=300)

#Table_markers
Idents(HIF) <- 'fig1'
Signatures_HIF<-data.frame()

for (i in 1:length(levels(HIF))) 
  {
z <- FindConservedMarkers(HIF,ident.1 = sort(levels(HIF))[i], grouping.var="orig.ident", logfc.threshold = 0.25, min.pct = 0.5, verbose = T, only.pos = T, assay = 'integrated')
z$gene <-rownames(z)
z$NCBI_id <- mapIds(org.Mm.eg.db, z$gene, 'ENTREZID', 'SYMBOL')
z$cluster <- sort(levels(HIF))[[i]]
Signatures_HIF<- rbind(Signatures_HIF, z)
}
Supplementary_Table_4 <- Signatures_HIF
write.csv(x=Signatures_HIF, file = "Supplementary_Table_4.csv", row.names = F) 

#Gene deregulated 
Idents(HIF)<-"fig1"
PTENvsHIF_genes <-data.frame()
HIF$celltype.treat <- paste(Idents(HIF), HIF$orig.ident, sep = "_")
Idents(HIF) <- "celltype.treat"
new.cluster.ids <- sort(c(levels(HIF@active.ident)[1:10]))
new.cluster.ids2 <- sort(c(levels(HIF@active.ident)[11:20]))

for (i in 1:length(new.cluster.ids)) 
  {
z = FindMarkers(HIF, ident.1 = new.cluster.ids[i], ident.2 = new.cluster.ids2[i], verbose = TRUE, logfc.threshold = 0.2, assay='RNA')
  z$gene <-rownames(z)
  z$NCBI_id <- mapIds(org.Mm.eg.db, z$gene, 'ENTREZID', 'SYMBOL')
  z$cluster <- sort(levels(HIF$fig1))[i]
  PTENvsHIF_genes<- rbind.fill(PTENvsHIF_genes, z)
  }
 write.csv(x=PTENvsHIF_genes, file = "HIFvsPTEN_DE_genes.csv", row.names = T) 

HIFvsPTEN_DE_genes <- PTENvsHIF_genes  
Idents(HIF)<-"fig1"
HIF_DE_genes_Table <- HIFvsPTEN_DE_genes[HIFvsPTEN_DE_genes$p_val_adj < 0.05 & HIFvsPTEN_DE_genes$avg_log2FC > 0.2, ]
x<- data.frame(table(HIF_DE_genes_Table$cluster))

x$Var1<-factor(x$Var1, levels=  levels(HIF$fig1))
my_color_palette <- hue_pal()(length(levels(HIF$fig1)))
ggplot(data = x, aes(x=Var1, y=Freq, fill=Var1))+geom_bar(stat="identity", color="black") + scale_fill_manual(values = my_color_palette) + theme(axis.title = element_text(size = 10), axis.line.x = element_line(size = 0), axis.line.y.left = element_line(size = 0.2), axis.text.y.left =  element_text(size = 8), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_line(size=0.2), axis.title.y = element_text(size = 10),axis.text.x = element_text(size = 8, angle=90, hjust = 1, vjust = 0.4, margin = margin(t =-0.1, unit =  "cm")), title = element_text(size = 8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + NoLegend() +ylab("number of\n downregulated genes")&ylim(0,300)
ggsave("Fig3B.eps", units = "mm", height = 50, width = 50, dpi=300)``

Supplementary_Table_5 <- PTENvsHIF_genes %>% filter(avg_log2FC > 0.2) %>% filter(p_val_adj < 0.05) %>% filter(cluster=="Luminal-C")
write.csv(x=Supplementary_Table_5 , file = "Supplementary_Table_5.csv", row.names = F) 

#In Luminal-C
Idents(HIF)<-"fig1"
z = FindMarkers(HIF, ident.1 = '3mo', ident.2 = '3moHIF', group.by='orig.ident', subset.ident="Luminal-C", logfc.threshold = 0.2, verbose = TRUE, assay='RNA', only.pos = T)
  z$gene <-rownames(z)
  z$NCBI_id <- mapIds(org.Mm.eg.db, z$gene, 'ENTREZID', 'SYMBOL')

  Gene_up <- z[z$p_val<0.05,]   
  KEGG_luminal_HIF <- enrichKEGG(Gene_up$NCBI_id, organism = "mmu", keyType="kegg", pvalueCutoff=0.05, pAdjustMethod = "fdr", universe,  minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.05)
KEGG_luminal_HIF<-setReadable(KEGG_luminal_HIF, org.Mm.eg.db, keyType = "ENTREZID")
write.table(KEGG_luminal_HIF@result, 'KEGG_HIF_Luminal.txt', sep = ' ')
my_gene_list <- KEGG_luminal_HIF@result[KEGG_luminal_HIF@result$Description=='HIF-1 signaling pathway',]
my_gene_list <- str_split(my_gene_list$geneID, "/", simplify=T)
Supplementary_Table_6 <- KEGG_luminal_HIF@result
write.csv(x = KEGG_luminal_HIF@result, file = "Supplementary_Table_6.csv", row.names = F)

#Figure 3
for (i in 1:(length(my_gene_list))) 
  {
VlnPlot(object = HIF, features = my_gene_list[[1,i]],pt.size = 0.1, assay = "RNA", idents=c("Luminal-C"), split.by = "orig.ident", cols = c('#F2756D','#21868F')) + theme(axis.title = element_text(size = 8), axis.line.x = element_line(size = 0), axis.line.y.left = element_line(size = 0.2), axis.text.y.left =  element_text(size = 8), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_line(size=0.2), axis.title.y = element_blank(), axis.text.x = element_text(size = 8, angle=0, hjust=0.5, vjust=0), title = element_text(size = 8)) + NoLegend()
name <- paste0(my_gene_list[[i]],"_HIF_Fig3.eps")
ggsave(name, units = "mm", height = 40, width = 40, dpi=300)  
 }

#Figure 3G
Figure3G <- c('Ifng', 'Tnf')
for (i in 1:(length(Figure3G))) 
  {
VlnPlot(HIF, features = Figure3G[[i]],pt.size = 0.1, assay = "RNA", split.by = "orig.ident", idents='T-cells', cols = c('#F2756D','#21868F')) + theme(axis.title = element_text(size = 8), axis.line.x = element_line(size = 0), axis.line.y.left = element_line(size = 0.2), axis.text.y.left =  element_text(size = 8), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_line(size=0.2), axis.title.y = element_blank(), axis.text.x = element_text(size = 8, angle=0, hjust=0.5, vjust=0), title = element_text(size = 8)) + NoLegend()
name <- paste0(Figure3G[[i]],"_HIF_Fig4D.eps")
ggsave(name, units = "mm", height = 40, width = 40, dpi=300)  
}

#Analysis of the Trajectory
#Seurat Object
Ech1 = CreateSeuratObject(data_10x,
                          project = "WT",
                          assay = "RNA",
                          min.cells = 10,
                          min.features = 100,
                          names.field = 1,
                          names.delim = "_",
                          meta.data = NULL
)

Ech1$orig.ident="1_WT"
Ech1[["percent.mt"]] <- PercentageFeatureSet(Ech1, pattern = "^mt-")
VlnPlot(Ech1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=1, ncol = 3)
Traj <- Ech1
Traj <- subset(Traj, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 50)
Traj <- NormalizeData(Traj)
Traj <- FindVariableFeatures(Traj, selection.method = "vst", nfeatures = 2000)

#Run the standard workflow for visualization and clustering
all.genes <- rownames(Traj)
Traj <- ScaleData(Traj, features = all.genes)
Traj <- RunPCA(Traj, npcs = 50, verbose = T)

# t-SNE and Clustering

Traj <- FindNeighbors(Traj, reduction = "pca")
Traj <- FindClusters(Traj,resolution =  0.6)
Traj <- RunTSNE(Traj, dims.use = 1:30, reduction.use = "pca", perplexity = 30)
Traj <- RunUMAP(Traj, dims = 1:30)

Traj@meta.data$Barcode <- rownames(Traj@meta.data)
Traj@meta.data <- right_join(clustres, Traj@meta.data, by='Barcode')
rownames(Traj@meta.data) <- Traj@meta.data$Barcode

Traj@meta.data <-Traj@meta.data %>% mutate(test=ifelse(grepl('-1', Barcode), '9mo', 
                                       ifelse(grepl('-2', Barcode), 'w15mo', 
                                              ifelse(grepl('-3', Barcode), '6mo', 
                                                     ifelse(grepl('-4', Barcode), '3mo', 'none')))))

Idents(Traj) <- 'Cluster'
Traj_subset <- subset(Traj, idents=c(10,12,20))
DimPlot(Traj_subset, split.by = 'test', order=c('3mo', '6', '9mo', '15'))

Idents(Traj_subset)<-'Cluster'
Idents(Traj_subset)<- factor(x = Idents(Traj_subset), levels = c(20,10,12))
z = FindMarkers(Traj_subset, ident.1 = "12", ident.2 = c("20", "10"), verbose = TRUE, logfc.threshold = 0.2, assay='RNA', only.pos = T)
  z$gene <-rownames(z)
  z$NCBI_id <- mapIds(org.Mm.eg.db, z$gene, 'ENTREZID', 'SYMBOL')
write.csv(x = z, file = "Supplementary_Table_.csv", row.names = F)
remove(z, Gene_up)

DimPlot(Traj_subset, label = F,label.size = 3, cols = c('#C0B6DB', '#78D7C8', '#FBFF6E'), order=c('3mo', '6', '9mo', '15', order = T), pt.size = 0.3, split.by = 'test')&xlim(c(-3, 12)) &ylim(c(-12, 5)) & theme(axis.text.x = element_text(size = 8, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8), title = element_text( size = 7), legend.title = element_text(size = 8), legend.key.height = unit(0.4, "lines"), legend.key.width = unit(0.5, "lines"), legend.text = element_text(size = 8))
ggsave("Fig6A.eps", units = "mm", height = 50, width = 100, dpi=300)

  FeaturePlot(Traj_subset, 'Hif1a', keep.scale = 'all', slot = 'data', label = F,label.size = 3, pt.size = 0.2, split.by = 'test', cols = c('gray', 'red'),order = T)&xlim(c(-3, 12)) &ylim(c(-12, 5))&theme(legend.position = ('bottom'))& theme(axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 1), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8), title = element_text( size = 7), legend.title = element_text(size = 8), legend.key.height = unit(0.4, "lines"), legend.key.width = unit(0.5, "lines"), legend.text = element_text(size = 8))
ggsave("Fig6D.eps", units = "mm", height = 50, width = 100, dpi=300)

my_gene_list <- c('Hif1a', 'Slc2a1', 'Hk2', 'Car9', 'Vegfa', 'Tgm2')
for (i in 1:(length(my_gene_list))) {
  
VlnPlot(Traj_subset, features = my_gene_list[i], cols = c('#C0B6DB', '#FBFF6E', '#78D7C8'), pt.size = 0.05)& theme(axis.title = element_text(size = 8), title = element_text(  size = 7), axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 1, angle = 0), axis.text.y = element_text(size = 8), legend.title = element_text(size = 8), legend.key.height = unit(0.4, "lines"), legend.key.width = unit(0.5, "lines"), legend.text = element_text(size = 8))     
name <- paste0(my_gene_list[[i]],"_traj.eps")
ggsave(name, units = "mm", height = 40, width = 45, dpi=300)  
}

#correlation TCGA
#a gene of interest

a <- log2(a)
b <- log2(as.numeric(graph$HIF1A))
cor<-cor.test(a, b, method = 'spearman')
cor2 <- cor.test(a, b, method = 'pearson')
plot <- data.frame(a, b)
ggplot(plot, (aes(x=as.numeric(a), y=as.numeric(b)))) +geom_point(stat="identity",shape=1, size=1, colour='black') + geom_smooth(method = "lm", fill="blue", colour="blue")+ theme(axis.text.x = element_text(size = 8, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8), title = element_text( size = 7), legend.title = element_text(size = 8), legend.key.height = unit(0.4, "lines"), legend.key.width = unit(0.5, "lines"), legend.text = element_text(size = 8), axis.line.x = element_line(colour = 'black', size = 0.5),axis.line.y = element_line(colour = 'black', size = 0.5), panel.background = element_rect(fill='white', colour='black', size=0, linetype = 1, color = 'black')) + xlab(label = data[i]) + ylab(label = 'TGM2')+ggtitle(cor[["p.value"]])
pval <- c(pval, cor[["p.value"]])
rho1 <- c(rho1, cor[["estimate"]])
pval2 <- c(pval2, cor2[["p.value"]])
rho2 <- c(rho2, cor2[["estimate"]])
gene <- c(gene, colnames(test)[i])
name <- paste0("HiF_finalFigure/old/correlation/", colnames(graph)[[i]],"_PTEN.eps")
ggsave(filename = name, units = "mm", height = 45, width = 45, dpi=300)  
}

#################################
#Trajectory inference using dyno
#Author : Celine Keime
#Date : 2020/08/24
#################################



#load libraries
###############

library(dyno)
library(tidyverse)
library(dplyr)
library(Seurat)



#parameters 
#############

#path to directory containing barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files
data.dir="S19035_S19177_S19192/outs/filtered_feature_bc_matrix"

#detection threshold (ie min UMI number to consider a gene is expressed)
detection = 5

#include genes detected in at least this number of cells
min.cells = 10

#include cells where at least this number of genes are detected
min.features = 100

#include cells with unique gene count less than this number
max.features=5000

#include cells with less than this percentage of mitochondrial reads
max.percent.mt = 20

#selected barcodes
clustres = read.table("S19035_S19177_S19192/outs/analysis/clustering/graphclust/clusters.csv",sep=",", header=TRUE, stringsAsFactors = FALSE)
clustres.selectedbarcodes = clustres[clustres$Cluster==10 | clustres$Cluster==12 | clustres$Cluster==20,  ]

#grouping labels
grouping_labels = c("3_9m","4_15m","2_6m","1_3m")
names(grouping_labels) = as.character(1:4)

#trajectory method
trajectory_method="slingshot"

#milestone id that is at the start of the trajectory
root_milestone_id = "3"

#size of dots representing cells
size_cells = 1



#data preprocessing
###################

#read 10X data
data_10x = Read10X(data.dir = data.dir)

#create Seurat object
data = CreateSeuratObject(counts=data_10x, min.cells=min.cells, min.features=min.features)

#percentage mitochondrial reads
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")

#filter data to exclude cells with more than max.percent.mt mitochondrial reads
#and to keep a number of features between min.features and max.features
data = subset(data, subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mt < max.percent.mt)

#normalize data
data = NormalizeData(data, normalization.method = "LogNormalize")

#filter data to exclude genes detected (with UMI > detection) in less than min.cells cells
keptgenes = which(rowSums( as.matrix(data@assays$RNA@counts) >= detection) >=  min.cells )
data = subset(data, features=keptgenes)

#retain only selected cells
selectedbarcodes = clustres.selectedbarcodes$Barcode
data = subset(data, cells=selectedbarcodes)



#data preparation for dyno
##########################

#genes in column
data_counts = Matrix::t(as(as.matrix(data@assays$RNA@counts), 'sparseMatrix'))/log(2)
data_expression = Matrix::t(as(as.matrix(data@assays$RNA@data), 'sparseMatrix'))/log(2)

dataset = wrap_expression(
  counts = data_counts, #counts values of genes (columns) within cells (rows)
  expression = data_expression #normalized expression values of genes (columns) within cells (rows)
)

#add grouping information according to sample number 
data_grouping = substr(rownames(data_counts),18,18) #sample number 
#use grouping_labels to replace sample number
for (l in unique(data_grouping)){
  data_grouping[data_grouping==l] = grouping_labels[l]
}
dataset = add_grouping(
  dataset,
  grouping = data_grouping
)

  
  
#trajectory inference
#####################

set.seed(1)
model = infer_trajectory(dataset, method=trajectory_method)
#rooting
model = model %>% add_root(root_milestone_id = root_milestone_id)



#plot trajectory 
#################

#colored by samples
plot_dimred(model, grouping = dataset$grouping, size_cells = size_cells)


#colored by clusters
data_grouping_cluster=NULL #grouping information according to cluster number
for (b in rownames(data_counts)){
  data_grouping_cluster = append(data_grouping_cluster, as.character(clustres.selectedbarcodes$Cluster[clustres.selectedbarcodes$Barcode==b]))
}
names(data_grouping_cluster) = rownames(data_counts)

plot_dimred(model, grouping = data_grouping_cluster, size_cells = size_cells)


  
  
#calculate feature importance for the whole trajectory
######################################################

overall_feature_importances = dynfeature::calculate_overall_feature_importance(model, expression_source = dataset$expression)
  
#file containing feature importance
write.table(overall_feature_importances, file="overall_feature_importances.txt", quote=FALSE, sep="\t", row.names=FALSE)

