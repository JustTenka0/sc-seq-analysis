library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

################################################################################################################
# 1. Upload Dataset
################################################################################################################
load("SRA713577_SRS3363004.sparse.RData")

# Initial filtering and adjusting rownames
sm@Dimnames[[1]] <- sub("_[^_]*$", "", sm@Dimnames[[1]])
genenames <- rownames(sm)
rownames(sm) <- make.names(genenames, unique=TRUE)

# pbmc is my dataset, Peripheral Blood Mesenchimal Cells
pbmc <- CreateSeuratObject(counts = sm, project = "PBMC_SRA713577",
                                         min.cells = 3, min.features = 200)

pbmc # 3568 samples
################################################################################################################
# 2. Quality control metrics
################################################################################################################
# For unknown reasons,  MT-proteins are indicated with MT. ("^MT\\.") and not MT- ("^MT-") 
# we can still change with rownames(sm) <- make.names(genenames, unique=TRUE)

# 1. Check with grep of MT e RBP if there are these genes
# grep("^MT-",rownames(pbmc),value = TRUE)
# Mitochondrial
grep("^MT\\.", rownames(pbmc), value = TRUE)
# grep("^MT-", rownames(pbmc), value = TRUE)   

# Ribosomal
grep("^RP[LS]", rownames(pbmc), value = TRUE)

# 2. Precentages computation
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT\\.")
pbmc[["percent.rbp"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[LS]")

# 3. Check if metadata are realistic (e.g. between 0 and 10-15 for mt)
head(pbmc@meta.data)
head(pbmc_before_filtering@meta.data)


# Do this command just one time to have a before filtering dataset
pbmc_before_filtering <- pbmc

# 4. Violin Plot (before filtering)
qcvlnplot_before <- VlnPlot(pbmc_before_filtering, ncol = 4, pt.size = 0, cols = "67", layer = "counts",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"))  & theme(
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),          
          axis.title.x = element_blank(),                                  
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 10, face = "italic")
        )

qcvlnplot_before <- qcvlnplot_before + plot_annotation(
          title = "Quality metrics before Filtering",
          theme = theme(
            plot.title = element_text(face = "bold", size = 16, hjust = 0.5, vjust = 1),
            plot.subtitle = element_text(size = 11, hjust = 0.5, face = "italic")
          )
        )

qcvlnplot_before

# Checking features
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = "67")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "67")
plot1 + plot2

plot3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.rbp", cols = "67")
plot3

################################################################################################################
# 3. Quality control & Filtering
################################################################################################################
# 1. Choosing thresholds
# Classic hard threshold (in Seurat tutorial)
lower_limit_feat = 200
upper_limit_feat = 2500
upper_limit_mt = 5


### Alternative: threshold computation, not needed ###
# stats_features <- pbmc@meta.data$nFeature_RNA
# median_feat <- median(stats_features)
# mad_feat    <- mad(stats_features)

## nFeatureRNA limits ##
# lower_limit_feat <- median_feat - (3 * mad_feat)
# upper_limit_feat <- median_feat + (3 * mad_feat)
# lower_limit_feat <- max(lower_limit_feat, 200)


## mt % limit ##
# stats_mt <- pbmc@meta.data$percent.mt
# median_mt <- median(stats_mt)
# mad_mt    <- mad(stats_mt)
# upper_limit_mt <- median_mt + (3 * mad_mt)



cat("----- Thresholds--------\n")
cat("nFeatureRNA (> 200): ",lower_limit_feat,"\n")
cat("nFeatureRNA (< 2500): ",upper_limit_feat,"\n")
cat("mt % (< 5): ",upper_limit_mt,"\n")

# 2. Apply thresholds
pbmc <- subset(pbmc, subset = nFeature_RNA > lower_limit_feat & 
                 nFeature_RNA < upper_limit_feat & 
                 percent.mt < upper_limit_mt)

pbmc

# 3. Update on the Violin plot
qcvlnplot_after <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), 
        ncol = 4, pt.size = 0, cols = "67") & theme(
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5), 
          axis.title.x = element_blank(),                                   
          axis.text.x = element_blank(),
          axis.title.y = element_text(size = 10, face = "italic")
        )

qcvlnplot_after <- qcvlnplot_after + plot_annotation(
  title = "Quality metrics after Filtering",
  theme = theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, vjust = 1),
    plot.subtitle = element_text(size = 11, hjust = 0.5, face = "italic")
  )
)

qcvlnplot_after

dim(pbmc_before_filtering)
dim(pbmc)


################################################################################################################
# 4. Normalization
################################################################################################################
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc@assays
pbmc@assays$RNA

################################################################################################################
# 5. Gene expression exploration
################################################################################################################
apply(pbmc[["RNA"]]$data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression

head(gene.expression, n=50)
VlnPlot(pbmc, features = c("MALAT1","GAPDH"), pt.size = 0, cols = "67")


cc.genes.updated.2019
CellCycleScoring(pbmc, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> pbmc

pbmc[[]]

#the default method -vst- computes (or better, estimates) the mean-variance relationship of each gene, and chooses the 2000 genes with hte highest variance. 
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
top10

all.genes <- rownames(pbmc)


pbmc <- ScaleData(pbmc, features = all.genes)
pbmc@assays$RNA

################################################################################################################
# 6. PCA
################################################################################################################
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# 1. Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:2, nfeatures = 10)
print(pbmc[["pca"]], dims = 1:10, nfeatures = 5)

# 2. Visualize first 2 PCs
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# 3. Visualize Elbow plot for choosing the right PCs
ElbowPlot(pbmc, ndims=35) +
  labs(
    title = "Elbow Plot - PCA Selection",
    x = "Principal Components (PCs)",
    y = "Standard Deviation"
  ) 

# 4. Choosing the right PC through the formula 
pc.touse <- (pbmc$pca@stdev)^2
pc.touse <- pc.touse/sum(pc.touse)
pc.touse <- cumsum(pc.touse)[1:50]
pc.touse <- min(which(pc.touse>=0.75))
pc.touse # 20

## If needed, re-calculate the percentage of variance for each PC (e.g. I have 10 PCs, I'll find a 65% variance)
pc.var <- (pbmc$pca@stdev)^2
pc.per <- pc.var / sum(pc.var)
pc.cum <- cumsum(pc.per)
variance_percentage <- pc.cum[10] * 100
print(variance_percentage)

################################################################################################################
# 7. Clustering
################################################################################################################
# 1. Configuration creation (choose how you feel the PCs and resolution0)
# === CONF. A (10 PC, Res 0.5) ===
pbmc10 <- FindNeighbors(pbmc, dims = 1:10)
pbmc10_05 <- FindClusters(pbmc10, resolution = 0.5)
pbmc10_05 <- RunUMAP(pbmc10_05, dims = 1:10)

table(Idents(pbmc10_05)) # to count the size of each cluster


# === CONF. B (15 PC, Res 0.8) ===
pbmc15 <- FindNeighbors(pbmc, dims = 1:15)
pbmc15_08 <- FindClusters(pbmc15, resolution = 0.8)
pbmc15_08 <- RunUMAP(pbmc15_08, dims = 1:15)

# === CONF. C (20 PC, Res 0.5) ===
pbmc20 <- FindNeighbors(pbmc, dims = 1:20)
pbmc20_05 <- FindClusters(pbmc20, resolution = 0.5)
pbmc20_05 <- RunUMAP(pbmc20_05, dims = 1:20)

table(Idents(pbmc20_05)) # to count the size of each cluster


# === CONF. D (20 PC, Res 0.8) ===
pbmc20 <- FindNeighbors(pbmc, dims = 1:20)
pbmc20_08 <- FindClusters(pbmc20, resolution = 0.8)
pbmc20_08 <- RunUMAP(pbmc20_08, dims = 1:20)
table(Idents(pbmc20_08))


# === CONF. E (30 PC, Res 0.8) (85% variance) === 
pbmc30 <- FindNeighbors(pbmc, dims = 1:30)
pbmc30_08 <- FindClusters(pbmc30, resolution = 0.8)
pbmc30_08 <- RunUMAP(pbmc30_08, dims = 1:30)


# 2. Plot creation
plot1 <- DimPlot(pbmc10_05, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "10 PCs - Res 0.5")
plot2 <- DimPlot(pbmc15_08, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "15 PCs - Res 0.8")
plot3 <- DimPlot(pbmc20_05, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "20 PCs - Res 0.5")
plot4 <- DimPlot(pbmc20_08, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "20 PCs - Res 0.8")
plot5 <- DimPlot(pbmc30_08, reduction = "umap", label = TRUE) + NoLegend() + labs(title = "30 PCs - Res 0.8") # 85%

# 3. Plot visualization, single or already paired
plot5
plot4

plot3+plot4
plot2+plot3

plot1 + plot3

plot1
plot3


# Optimal
head(pbmc20_05[[]],5)
DimPlot(pbmc20_05, reduction = "pca")


pbmc <- RunTSNE(pbmc, dims=1:10)
DimPlot(pbmc, reduction = "tsne")

pbmc <- RunUMAP(pbmc, dims = 1:10)
#if you cannot install UMAP, t_SNE is anyway ok for your project!
DimPlot(pbmc, reduction = "umap")

################################################################################################################
# 8. Find the Markers
################################################################################################################
# Find all markers of cluster 2 versus all the others 
# (helpful if you want to see top expressed in a cluster later on)
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25, test.use = "wilcox")
head(cluster2.markers, n = 5)

# Find top5 marker for each cluster in the chosen dataset 
pbmc.markers <- FindAllMarkers(pbmc20_08, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
filtered_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  ungroup()

# Explore markers
print(filtered_markers, n = 65)
View(filtered_markers)
# Extract the markers of a subset of a cluster
subset(filtered_markers, cluster == "10")$gene

################################################################################################################
# 9. Marker plots
################################################################################################################
# FEATURE PLOT with a bilanced marker genes choice
FeaturePlot(pbmc20_08, features = c("IL7R","CD8A", "KLRF1", "LYZ", "MS4A7", "LILRA4","PPBP",  "TCL1A", "IGHG2"), ncol = 3)

# DOT PLOT
DotPlot(pbmc20_08,
        features = c("CCR7","IL7R","CD8A", "GZMH", "GZMK", "CD8B", "KLRF1","LYZ", "MS4A7", "HBB" ,  "TCL1A", "IGHG2", "CLDN5","PPBP","FCGR3B", "LILRA4")) +
  theme(
          # Inclina il testo dell'asse X (i geni) di 45 gradi e lo allinea a destra (vjust/hjust)
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 11, face = "bold"),
          
          # Opzionale: Aumenta lo spazio/margine inferiore del plot per non far tagliare i nomi dei geni lunghi
          plot.margin = margin(t = 10, r = 10, b = 20, l = 10, unit = "pt")
        )




# HEATMAP
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc20_08, features = top10$gene) + NoLegend()

################################################################################################################
# 10. Cluster genes analysis
################################################################################################################
# cluster 3
# Find specific markers for cluster 3 using the Wilcoxon Rank Sum test
cluster3.markers <- FindMarkers(pbmc20_08, ident.1 = 3, min.pct = 0.25, test.use = "wilcox")
# Sort markers by average log2 fold-change in descending order (highest expression first)
cluster3.markers <- cluster3.markers[order(-cluster3.markers$avg_log2FC), ]
# Display top 10 most upregulated genes for Cluster 3
head(cluster3.markers, n = 10)


# cluster 10
# Find specific markers for cluster 10
cluster10.markers <- FindMarkers(pbmc20_08, ident.1 = 10, min.pct = 0.25, test.use = "wilcox")
# Sort markers by average log2 fold-change in descending order
cluster10.markers <- cluster10.markers[order(-cluster10.markers$avg_log2FC), ]
# Display the full dataframe of the top 15 markers (fixed from just printing rownames)
head(cluster10.markers, n = 15)


# Similar cluster 4,9
# 1. Isolate markers specifically enriched in Cluster 9 (Naive B Cells)
cluster9.markers <- FindMarkers(pbmc20_08, ident.1 = 9, min.pct = 0.25, test.use = "wilcox")
cluster9.markers <- cluster9.markers[order(-cluster9.markers$avg_log2FC), ]
head(cluster9.markers, n = 10) # Fixed: properly printing cluster 9 instead of cluster 6

# 2. Identify shared markers for both B-cell/Plasma related clusters (4 and 9) against all other cells
cluster4AND9.markers <- FindMarkers(pbmc20_08, ident.1 = c(4, 9), min.pct = 0.25, test.use = "wilcox")
cluster4AND9.markers <- cluster4AND9.markers[order(-cluster4AND9.markers$avg_log2FC), ]
head(cluster4AND9.markers, n = 10)

# 3. Direct pairwise Differential Expression (DE): Cluster 9 versus Cluster 4 (Plasma Cells)
# Positive log2FC indicates genes upregulated in cluster 9; negative log2FC indicates upregulation in cluster 4
cluster9vs4.markers <- FindMarkers(pbmc20_08, ident.1 = 9, ident.2 = 4, min.pct = 0.25, test.use = "wilcox")
cluster9vs4.markers <- cluster9vs4.markers[order(-cluster9vs4.markers$avg_log2FC), ]

# Display top 10 genes significantly upregulated in Cluster 9 (Positive log2FC)
head(cluster9vs4.markers, n = 10)

# Display top 10 genes significantly upregulated in Cluster 4 (Negative log2FC)
tail(cluster9vs4.markers, n = 10)

VlnPlot(pbmc20_05, features = c("IGHG2", "S100A4"), pt.size = 0)

################################################################################################################
# 11. Assign cluster names
################################################################################################################

# Biological names assigned
new.cluster.ids <- c(
  "Naive CD4+ T Cells",                # Cluster 0
  "Memory CD4+ T Cells",               # Cluster 1
  "NK Cells",                          # Cluster 2
  "Erythroid Cells",                   # Cluster 3 (Contaminazione HBB)
  "Plasma Cells",                      # Cluster 4 (Fabbriche di IGHG)
  "Cytotoxic T Cells (GZMH+)",         # Cluster 5 ()
  "Cytotoxic T Cells (GZMK+)",      # Cluster 6 ()
  "Cytotoxic T Cells (CD8B+)",         # Cluster 7 (CD8B+)
  "Monocytes",                   # Cluster 8 (LYZ+/MAFB+)
  "Naive B Cells",                     # Cluster 9 (TCL1A+/IGHD+)
  "Endothelial Cells /\
   Platelet Aggregates",                 # Cluster 10 (CLDN5+)
  "Neutrophils",                       # Cluster 11 (FCGR3B+)
  "Plasmacytoid Dendritic Cells"       # Cluster 12 (SCT+)
)


names(new.cluster.ids) <- levels(pbmc20_08)
pbmc20_08 <- RenameIdents(pbmc20_08, new.cluster.ids)
names(new.cluster.ids) <- levels(pbmc20_08)
pbmc20_08 <- RenameIdents(pbmc20_08, new.cluster.ids)

# Final UMAP
DimPlot(pbmc20_08, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 4) + NoLegend()


