#https://satijalab.org/seurat/articles/get_started.html
#
library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
#The . values in the matrix represent 0s (no molecules detected). Since most values in an scRNA-seq matrix are 0, Seurat uses a sparse-matrix representation whenever possible. This results in significant memory and speed savings for Drop-seq/inDrop/10x data.}
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#The number of unique genes and total molecules are automatically calculated during CreateSeuratObject()
#You can find them stored in the object meta data
head(pbmc@meta.data, 5)

#Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in
#pbmc[["RNA"]]@data.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#es lo mismo que:   pbmc <- NormalizeData(pbmc)

#Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10
# plot variable features 
#en dataset pequeño
data("pbmc_small")
pbmc_small
class(pbmc_small)
VariableFeaturePlot(object = pbmc_small)

#en este ejemplo
VariableFeaturePlot(object = pbmc)

#Scaling the data for dimension reduction
#results of this are stored in pbmc[["RNA"]]@scale.data
#solo en los 2000 genes seleccionados arriba

pbmc <- ScaleData(pbmc)
#en todos los genes:
#all.genes <- rownames(pbmc)
#pbmc <- ScaleData(pbmc, features = all.genes)

#Perform linear dimensional reduction
#PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:6, cells = 500, balanced = TRUE)

#‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line). In this case it appears that there is a sharp drop-off in significance after the first 10-12 PCs.
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

#An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot() function). In this example, we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.
ElbowPlot(pbmc)

#Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
#we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the FindNeighbors() function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs).

#To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default) or SLM [SLM, Blondel et al., Journal of Statistical Mechanics], to iteratively group cells together, with the goal of optimizing the standard modularity function. The FindClusters() function implements this procedure, and contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. 
head(Idents(pbmc), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
#pip install umap-learn). Details on this package can be found here: https://github.com/lmcinnes/umap. 
pbmc <- RunUMAP(pbmc, dims = 1:10)

#Guardar hasta acá
saveRDS(pbmc, file = "pbmc_tutorial.rds")

#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
#pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#pbmc.markers %>%
 # group_by(cluster) %>%
 # slice_max(n = 2, order_by = avg_log2FC)

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("MS4A1", "CD79A"), slot = "counts", log = TRUE)


#Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet", "??")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "pbmc_tutorial.rds")

##############
#Otros datos en librería SeuratData
#https://github.com/satijalab/seurat-data
library(SeuratData)
InstallData("pbmc3k")
data("pbmc3k")

#################
#Datos simulados
#### Librerias ####
# RNA-seq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scater")
BiocManager::install("splatter")
BiocManager::install("cowplot")

library("scater")
library("splatter")
# Plotting
library("cowplot")
# Tidyverse
library("tidyverse")#install.packages("tidyverse")

#### Simulacion de población con diferentes grupos de celulas ####
sim.groups <- splatSimulateGroups(batchCells = 1000,
                                  group.prob = c(0.2, 0.2, 0.2, 0.2, 0.2),                      de.prob = c(0.05, 0.05, 0.08, 0.05, 0.06),                      de.facLoc = 0.1,
                                  de.facScale = 0.4,
                                  dropout.type = "experiment",
                                  seed = 1)
sim.groups <- logNormCounts(sim.groups)

#### Matriz de conteos ####
matriz_conteos <- sim.groups@assays@data@listData[["logcounts"]]
write.csv(matriz_conteos, file = "matriz_conteos.csv")
dim(matriz_conteos)
#write.csv(matriz_conteos, file = "matriz_conteos.csv")
?CreateSeuratObject
pbmc <- CreateSeuratObject(counts = matriz_conteos)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in
#pbmc[["RNA"]]@data.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#es lo mismo que:   pbmc <- NormalizeData(pbmc)

#Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the most highly variable genes
top20 <- head(VariableFeatures(pbmc), 20)
top20
#en este ejemplo
VariableFeaturePlot(object = pbmc)

#Scaling the data for dimension reduction
#results of this are stored in pbmc[["RNA"]]@scale.data
#solo en los 2000 genes seleccionados arriba

pbmc <- ScaleData(pbmc)
#en todos los genes:
#all.genes <- rownames(pbmc)
#pbmc <- ScaleData(pbmc, features = all.genes)

#Perform linear dimensional reduction
#PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 100, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:6, cells = 100, balanced = TRUE)

#‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line). In this case it appears that there is a sharp drop-off in significance after the first 10-12 PCs.
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

#An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot() function). In this example, we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.
ElbowPlot(pbmc)
#5 PC es suficiente en este ejemplo

#Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:5)
pbmc <- FindClusters(pbmc, resolution = 0.5)
#we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the FindNeighbors() function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs).

#To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default) or SLM [SLM, Blondel et al., Journal of Statistical Mechanics], to iteratively group cells together, with the goal of optimizing the standard modularity function. The FindClusters() function implements this procedure, and contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. 
head(Idents(pbmc), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
#pip install umap-learn). Details on this package can be found here: https://github.com/lmcinnes/umap. 
pbmc <- RunUMAP(pbmc, dims = 1:10)

#Guardar hasta acá
saveRDS(pbmc, file = "pbmc_tutorial.rds")

#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

VlnPlot(pbmc, features = c("Gene464", "Gene735"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("Gene464", "Gene735"), slot = "counts", log = TRUE)



DimPlot(pbmc, reduction = "umap", pt.size = 0.5) + NoLegend()

saveRDS(pbmc, file = "pbmc_tutorial.rds")


