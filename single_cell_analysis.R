

#######################################################################
# Setting up packages                                                 #
#######################################################################
#                                                                     #

# set up the packages, use the "if(cond) expr" construct to check if the package has already been installed to avoid re-stalling them.
# load the packages that we need to use.


if(!require(data.table)) install.packages("data.table")
if(!require(Seurat)) install.packages("Seurat")
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(Matrix)) install.packages("Matrix")
if(!require(devtools)) install.packages("devtools")
if(!require(BiocManager)) install.packages("BiocManager")
if(!require(dplyr)) install.packages("dplyr")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(kableExtra)) install.packages("kableExtra")

library(data.table)
library(Seurat)
library(tidyverse)
library(Matrix)
library(dplyr)
library(ggplot2)
library(kableExtra)

sessionInfo()

#                                                                     #
#                                                                     #
#######################################################################





#######################################################################
# Downloading datasets                                                #
#######################################################################
#                                                                     #

# download datasets
# find information of GSE10086 from the webpage "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866"
# download 3 csv datasets we need into local working direction (in the format of .gz compressed files)
# read the data from the downloaded files and save as .Rdata objects
# load the .Rdata objects
# use the if...else... control-flow construct to check if the .Rdata files already exist before re-downloading them.

if(file.exists("CBMC8K_RNA_umi.Rdata")) load("CBMC8K_RNA_umi.Rdata") else{
  url_CBMC8K_RNA_umi <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DRNA%5Fumi%2Ecsv%2Egz"
  download.file(url_CBMC8K_RNA_umi, destfile = "./CBMC8K_RNA_umi.csv.gz")
  CBMC8K_RNA_umi <- read.csv("CBMC8K_RNA_umi.csv.gz", header=T, row.names=1)
  save(CBMC8K_RNA_umi, file = "CBMC8K_RNA_umi.Rdata")
  load("CBMC8K_RNA_umi.Rdata")
}

if(file.exists("CBMC8K_ADT_umi.Rdata")) load("CBMC8K_ADT_umi.Rdata") else{
  url_CBMC8K_ADT_umi <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100866&format=file&file=GSE100866%5FCBMC%5F8K%5F13AB%5F10X%2DADT%5Fumi%2Ecsv%2Egz"
  download.file(url_CBMC8K_ADT_umi, destfile = "./CBMC8K_ADT_umi.csv.gz")
  CBMC8K_ADT_umi <- read.csv("CBMC8K_ADT_umi.csv.gz", header=T, row.names=1)
  save(CBMC8K_ADT_umi, file = "CBMC8K_ADT_umi.Rdata")
  load("CBMC8K_ADT_umi.Rdata")
}


dim(CBMC8K_RNA_umi) # 36280  8617
dim(CBMC8K_ADT_umi) # 13 8617

#                                                                     #
#                                                                     #
#######################################################################





#######################################################################
# Exploring the dataset: separating human and mouse cells             #
#######################################################################
#                                                                     #

length(rownames(CBMC8K_RNA_umi)) # 36280
n_distinct(rownames(CBMC8K_RNA_umi)) # 36280

rownames(CBMC8K_RNA_umi)[1:10]
rownames(CBMC8K_RNA_umi)[36271:36280]

species <- tidyr::separate(as.data.frame(rownames(CBMC8K_RNA_umi)),col="rownames(CBMC8K_RNA_umi)", into=c("species","gene"),sep="_")
dim(species) # 36280     2

sum(species$species=="HUMAN") # 20400
sum(species$species=="MOUSE") # 15879
#20400+15879 # = 36279
index_human <- species$species=="HUMAN"
index_mouse <- species$species=="MOUSE"
index_human_or_mouse <- species$species=="HUMAN"|species$species=="MOUSE"
sum(index_human_or_mouse) # 36279

CBMC8K_RNA_umi_human <- CBMC8K_RNA_umi[index_human,]
dim(CBMC8K_RNA_umi_human) # 20400  8617
CBMC8K_RNA_umi_human[1:10,1:10]
CBMC8K_RNA_umi_human[100:110,523:533]

CBMC8K_RNA_umi_mouse <- CBMC8K_RNA_umi[index_mouse,]
dim(CBMC8K_RNA_umi_mouse) # 15879  8617
CBMC8K_RNA_umi_mouse[1:10,1:10]

percent_human <- 0
for(i in 1:8617){
  percent_human[i] <- sum(CBMC8K_RNA_umi_human[,i])/(sum(CBMC8K_RNA_umi_human[,i])+sum(CBMC8K_RNA_umi_mouse[,i]))
}
percent_human
save(percent_human, file="percent_human.Rdata")

cell_human_index <- percent_human>0.9
cell_mouse_index <- percent_human<0.1
cell_mixed_index <- percent_human>0.1&percent_human<0.9
sum(cell_human_index) # 8005
sum(cell_mouse_index) # 579
sum(cell_mixed_index) # 33

CBMC8K_RNA_umi_human_cell_gene <- CBMC8K_RNA_umi[index_human,cell_human_index]
CBMC8K_ADT_umi_human <- CBMC8K_ADT_umi[,cell_human_index]

save(CBMC8K_RNA_umi_human_cell_gene, file = "CBMC8K_RNA_umi_human_cell_gene.Rdata")
save(CBMC8K_ADT_umi_human, file = "CBMC8K_ADT_umi_human.Rdata")

#                                                                     #
#                                                                     #
#######################################################################









#######################################################################
# Seurat package: human cells and genes                               #
#######################################################################
#                                                                     #

# Establish a Seurat object #
CBMC8K_human <- CreateSeuratObject(counts  = as.sparse(CBMC8K_RNA_umi_human_cell_gene),project = "CBMC8K")
CBMC8K_human # 20400 features across 8005 samples within 1 assay 

# to preview the data matrix:
GetAssayData(object = CBMC8K_human, slot = "counts")

# or:
CBMC8K_human[["RNA"]]@counts

# use @meta.data to have a brief look at the numbers of counts and features of each cell
CBMC8K_human@meta.data

# Visualize the distribution of feature numbers and total counts as violin plots
VlnPlot(object = CBMC8K_human, features = "nFeature_RNA", log = TRUE) # nGene
VlnPlot(object = CBMC8K_human, features = "nCount_RNA", log = TRUE) # nUMI

# Use FeatureScatter to visualize the relationship between nFeature and nCount 
FeatureScatter(object = CBMC8K_human, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")




# Quality control: filtering the single-cell RNA data #
# Cells with less than 500 genes detected and cells with a total number of UMIs or genes that is more than 3 s.d. above or below the mean (at log10 level) will be removed.

# convert nFeature and nCount to log10 level and add to @meta.data
CBMC8K_human[["nFeature_log10"]] <- log10(CBMC8K_human@meta.data$nFeature_RNA)
CBMC8K_human[["nCount_log10"]] <- log10(CBMC8K_human@meta.data$nCount_RNA)
CBMC8K_human@meta.data

summary(CBMC8K_human@meta.data$nFeature_log10)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.461   2.873   2.946   2.959   3.027   3.702 

summary(CBMC8K_human@meta.data$nCount_log10)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.870   3.198   3.312   3.320   3.415   4.492 

# calculate the upper and lower limit of nFeature
nFeature_upper <- mean(CBMC8K_human@meta.data$nFeature_log10) + 3*sd(CBMC8K_human@meta.data$nFeature_log10)
nFeature_lower <- mean(CBMC8K_human@meta.data$nFeature_log10) - 3*sd(CBMC8K_human@meta.data$nFeature_log10)
10^nFeature_upper # 2368.07
10^nFeature_lower # 349.5589
min_nFeature <- 500
max(min_nFeature,10^nFeature_lower) #500

# Visualize the limits of QC metrics as a violin plot
p_feature <- VlnPlot(object = CBMC8K_human, features = "nFeature_RNA", log = TRUE, do.return=TRUE) 
p_feature + 
  geom_hline(yintercept = 10^nFeature_upper, color="red")+
  geom_hline(yintercept = max(min_nFeature,10^nFeature_lower), color="red")

# do similar calculations for nCount
nCount_upper <- mean(CBMC8K_human@meta.data$nCount_log10) + 3*sd(CBMC8K_human@meta.data$nCount_log10)
nCount_lower <- mean(CBMC8K_human@meta.data$nCount_log10) - 3*sd(CBMC8K_human@meta.data$nCount_log10)
10^nCount_upper # 7617.333
10^nCount_lower # 572.6996

p_count <- VlnPlot(object = CBMC8K_human, features = "nCount_RNA", log = TRUE, do.return=TRUE) # nUMI
p_count + 
  geom_hline(yintercept = 10^nCount_upper, color="red")+
  geom_hline(yintercept = 10^nCount_lower, color="red")

# according to this we will filter cells and only keep the ones with nFeature of 500-2368 and nCount of 572-7617.
CBMC8K_human <- subset(CBMC8K_human, subset = nFeature_RNA > 500 & nFeature_RNA < 10^nFeature_upper & nCount_RNA > 10^nCount_lower & nCount_RNA < 10^nCount_upper)

CBMC8K_human # 20400 features across 7766 samples within 1 assay (239 cells were filtered out)
# 8005-7766 # =239

# visulize nFeature and nCount after filtering cells
VlnPlot(object = CBMC8K_human, features = "nFeature_RNA", log = TRUE)
VlnPlot(object = CBMC8K_human, features = "nCount_RNA", log = TRUE)
FeatureScatter(object = CBMC8K_human, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")



# Normalization #

# Now we will normalize the data. 
# By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in [["RNA"]]@data.
CBMC8K_human <- NormalizeData(CBMC8K_human)
# normalization.method = "LogNormalize" and scale.factor = 10000 are both default values

# these 2 lines of codes both give preview of normalized data 
CBMC8K_human[["RNA"]]@data
GetAssayData(object = CBMC8K_human)



# Identification of highly variable features (feature selection) #

# use "nfeatures" to define the number of features. here we use 1000.
CBMC8K_human <- FindVariableFeatures(CBMC8K_human, 
                                                selection.method = "vst", 
                                                nfeatures = 1000)

# have access to the list of highly variable features via "VariableFeatures(CBMC8K_human)"
length(VariableFeatures(CBMC8K_human))

# the 12 most highly variable genes
top12 <- head(VariableFeatures(CBMC8K_human), 12)
top12_name <- separate(as.data.frame(top12),1,into = c("discard","name"), sep = "-")[,2]

# plot variable features with and without labels
plot_1 <- VariableFeaturePlot(CBMC8K_human)
LabelPoints(plot = plot_1, points = top12, labels = top12_name, repel = TRUE)




# scaling #
# each gene is scaled to have mean expression of 0 and variance of 1.
CBMC8K_human <- ScaleData(CBMC8K_human, features = rownames(CBMC8K_human))

# preview of scaled data:
CBMC8K_human[["RNA"]]@scale.data[1:5, 1:5]

# or: 
GetAssayData(object = CBMC8K_human, slot = "scale.data")[1:5, 1:5]





# clustering #

# linear dimensional reduction using PCA #

# we will perform PCA on the scaled data.
# By default, only the previously determined highly variable features will be used as input.
CBMC8K_human <- RunPCA(CBMC8K_human, features = VariableFeatures(object = CBMC8K_human))

# by default RunPCA will extract the first 50 components (npcs = 50)
CBMC8K_human[["pca"]]
CBMC8K_human[["pca"]]@misc

# examine the first 5 components
print(CBMC8K_human[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize top genes associated with reduction components
plot_1 <- VizDimLoadings(CBMC8K_human, dims = 1, reduction = "pca",combine = TRUE)
plot_1 + theme(axis.text.x = element_text(angle = 45, hjust=1))
VizDimLoadings(CBMC8K_human, dims = 2, reduction = "pca",combine = TRUE)

# Graphs the output of a dimensional reduction technique on a 2D scatter plot
DimPlot(CBMC8K_human, reduction = "pca")

# In particular DimHeatmap allows for easy exploration of the primary sources of heterogeneity in a dataset, 
# and can be useful when trying to decide which PCs to include for further downstream analyses. 
# Both cells and features are ordered according to their PCA scores.
# Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum, which dramatically speeds plotting for large datasets.
DimHeatmap(CBMC8K_human, dims = 1, cells = 500, balanced = TRUE) 
DimHeatmap(CBMC8K_human, dims = 1:20, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset: JackStraw plot
CBMC8K_human <- JackStraw(CBMC8K_human, num.replicate = 100)
CBMC8K_human <- ScoreJackStraw(CBMC8K_human, dims = 1:20)
JackStrawPlot(CBMC8K_human, dims = 1:20)

# Determine the ‘dimensionality’ of the dataset: Elbow plot
ElbowPlot(CBMC8K_human)

# From the "elbow plot", the turning point is around PC11-17-ish. To be safe, I chose the first 17 principle components here for downstream clustering. 

# save the Seurat object for future use
saveRDS(CBMC8K_human, file = "CBMC8K_human.rds")

# When need to re-load the Seurat project, use the following code:
# CBMC8K_human <- readRDS("CBMC8K_human.rds")





# clustering #

# choose the first 17 PCs here
CBMC8K_human <- FindNeighbors(CBMC8K_human, dims = 1:17, k.param = 40)
CBMC8K_human <- FindClusters(CBMC8K_human, resolution = 0.8)
head(Idents(CBMC8K_human), 5)

# non-linear dimensional reduction (UMAP/tSNE)
# UMAP
CBMC8K_human <- RunUMAP(CBMC8K_human, dims = 1:17, min.dist = 0.75)
DimPlot(CBMC8K_human, reduction = "umap", label=T)

# TSNE
CBMC8K_human <- RunTSNE(seed.use = 1, object = CBMC8K_human,reduction = "pca",dims = 1:17)
DimPlot(CBMC8K_human, reduction = "tsne", label=T)

# visualize the clusters using PCA
DimPlot(CBMC8K_human, reduction = "pca", label=T)

saveRDS(CBMC8K_human, file = "CBMC8K_human.rds")




## Finding differentially expressed features (cluster biomarkers) ##

# find markers for every cluster compared to all remaining cells, report only the positive ones
CBMC8K_human.markers <- FindAllMarkers(CBMC8K_human, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# it takes some time to find all the markers so we'd better store the markers in .Rdata object in case we need to re-load them later.
save(CBMC8K_human.markers, file = "CBMC8K_human.markers.Rdata")
load("CBMC8K_human.markers.Rdata")

# examine the top 10 markers of each cluster
CBMC8K_human.markers %>% dplyr::filter(cluster==0) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==1) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==2) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==3) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==4) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==5) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==6) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==7) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==8) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==9) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==10) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==11) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==12) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==13) %>% top_n(n = 10, wt = avg_logFC)
CBMC8K_human.markers %>% dplyr::filter(cluster==14) %>% top_n(n = 10, wt = avg_logFC)

# Blood cells can be classified into different groups, characterized by combinations of common biomarkers.
# For example, CD8 T-cells have the feature of CD3+ CD4- CD8+, and B-cells are CD3- CD19+.
# We can visualize expression of marker genes in different clusters. 

# VlnPlot shows expression probability distributions across clusters
VlnPlot(CBMC8K_human, features = c("HUMAN-CD8A", "HUMAN-CD19"))

# you can plot raw counts as well
VlnPlot(CBMC8K_human, features = c("HUMAN-CD8A", "HUMAN-CD19"), slot = "counts", log = TRUE)

# FeaturePlot visualizes feature expression on a tSNE or PCA plot
FeaturePlot(CBMC8K_human, features = c("HUMAN-CD8A", "HUMAN-CD19"), reduction = "tsne")

# There are additional methods such as RidgePlot, and DotPlot.
RidgePlot(CBMC8K_human, features = c("HUMAN-CD8A", "HUMAN-CD19"))
plot_1 <- DotPlot(CBMC8K_human, features = c("HUMAN-CD8A", "HUMAN-CD19", "HUMAN-CD3E", "HUMAN-CD4", "HUMAN-CD34"))
plot_1 + theme(axis.text.x = element_text(angle = 45, hjust=1))

# DoHeatmap generates an expression heatmap for given cells and features. 

# Downsample the clusters to a maximum of 300 cells each (makes the heatmap easier to see for small clusters)
CBMC8K_human_small <- subset(CBMC8K_human, downsample = 300)

# We will use the top 10 markers (or all markers if less than 10) for each cluster to generate the heatmap.
top10 <- CBMC8K_human.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(CBMC8K_human_small, features = top10$gene)  + NoLegend()





# Assigning cell type identity to clusters #

# There are different types of blood cells, and we have accumulated significant amount of knowledge to recognize cell types based on biomarkers.
# Here's a table of canonical markers to easily match the unbiased clustering to known cell types:

celltype <- c("Naive CD4 T", "Memory CD4 T", "CD8 T", "B", "NK", "CD14+ Mono", "CD16+ Mono", "CD34+", "MK", "DC", "pDCs")
markergene <- c("CD3+ CD4+ CD8- CD2+/- CD57-", "CD3+ CD4+ CD8- CD2++ CD57-", "CD3+ CD4- CD8+", "CD3- CD19+", "CD16+ CD57+", "CD14++ CD16- CD15+", "CD14+ CD16+ CD15-", "CD34+", "CD41+", "CD141+ CD1c+", "CD11c+ CD45R+ CD123+ CD303+")
cell_type_table <- data.frame(CellType = celltype, Markers = markergene)
kable(cell_type_table) %>%
  kable_styling(bootstrap_options = "striped" , full_width = F , position = "center") %>%
  kable_styling(bootstrap_options = "bordered", full_width = F , position ="center") %>%
  column_spec(1,bold = F ) %>%
  column_spec(2,bold =T )


# Now we want to assign cell types to the clusters we got from the clustering process.
# We can visit the identities of each cell using "Idents(CBMC8K_human)", or "CBMC8K_human@meta.data$seurat_clusters". 
# levels(CBMC8K_human) returns the 15 clusters: "0" - "14"

Idents(CBMC8K_human)[1:10]
CBMC8K_human@meta.data$seurat_clusters[1:10]
levels(CBMC8K_human)

# Based on the markers of each cluster we identified earlier, and the knowledge of features of each cell type, we can now assign the cell types to these clusters.
new.cluster.ids <- c("Naive CD4 T","CD14+ Mono", "Memory CD4 T", "NK", "CD14+ Mono", "B", "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK", "Multiplets", "CD34+", "MK", "DC", "pDCs")
names(new.cluster.ids) <- levels(CBMC8K_human)
CBMC8K_human <- RenameIdents(CBMC8K_human, new.cluster.ids)

# Now we can see the IDs of each cell have been updated to their cell types.
Idents(CBMC8K_human) 
# "levels(CBMC8K_human)" returns the levels of these IDs (cell types).
levels(CBMC8K_human) 
# "CBMC8K_human@meta.data$seurat_clusters" still returns their original clusters.
CBMC8K_human@meta.data$seurat_clusters 

# If we want to reset the cell types, we can first use the following code to reset the IDs using the cluster numbers:
# Idents(CBMC8K_human) <- CBMC8K_human@meta.data$seurat_clusters

# Now we can add the cell type information to the tsne plot instead of the clustering numbers:
DimPlot(CBMC8K_human, reduction = "tsne", label = TRUE, pt.size = 0.5) 

# Always save the object so we can easily re-load it later.
saveRDS(CBMC8K_human, file = "CBMC8K_human.rds")
#                                                                     #
#                                                                     #
#######################################################################









#######################################################################
# Seurat package: Add protein data                                    #
#######################################################################
#                                                                     #

# Adding protein data to the Seurat object #

# Cell surface protein markers have also been accessed in this project. So we can intergrate this information to our Seurat project.
# By comparing to mouse cells as a non-specific background negative control, the authors ruled out 3 antibodies from their 13 antibodies used. So we first delete these 3 antibodies from the result.
load("CBMC8K_ADT_umi_human.Rdata")
CBMC8K_ADT_umi_human_10rows <- CBMC8K_ADT_umi_human[setdiff(rownames(x = CBMC8K_ADT_umi_human), c("CCR5", "CCR7", "CD10")), ]
dim(CBMC8K_ADT_umi_human_10rows) #  10 8005

# Because we filtered out some cells at the quality control stage, we now only keep the cells that exist in the final Seurat object.
CBMC8K_ADT_umi_human_10rows_7766 <- CBMC8K_ADT_umi_human_10rows[,colnames(CBMC8K_human)]
dim(CBMC8K_ADT_umi_human_10rows_7766) #  10 7766
identical(colnames(CBMC8K_human),colnames(CBMC8K_ADT_umi_human_10rows_7766))

# Now we add the protein data to the "ADT" assay of the Seurat object.
CBMC8K_human[["ADT"]] <- CreateAssayObject(counts = CBMC8K_ADT_umi_human_10rows_7766)
CBMC8K_human
saveRDS(CBMC8K_human, file = "CBMC8K_human.rds")

# Normalization: CLR
# For the protein data, we use the centered log ratio (CLR) transformation to normalize it.
CBMC8K_human <- NormalizeData(CBMC8K_human, assay = "ADT", normalization.method = "CLR")

# Scaling
CBMC8K_human <- ScaleData(CBMC8K_human, assay = "ADT")
saveRDS(CBMC8K_human, file = "CBMC8K_human.rds")

# We can check the original counts, normalized data and scaled data using the following codes:
CBMC8K_human[["ADT"]]@counts[1:5,1:5]
CBMC8K_human[["ADT"]]@data[1:5,1:5]
CBMC8K_human[["ADT"]]@scale.data[1:5,1:5]

# Now we have two assays of this Seurat object: a default "RNA" assay and a new "ADT" assay.
CBMC8K_human$RNA
CBMC8K_human$ADT

# We can use "DefaultAssay(CBMC8K_human)" to define which assay is considered as the default assay.
DefaultAssay(CBMC8K_human) <- "RNA"







# Visualize protein expression information #

# We can visualize the overall protein expression pattern in different clusters using "DotPlot".
plot_1 <- DotPlot(CBMC8K_human, assay = "ADT", features = c("adt_CD3", "adt_CD4", "adt_CD11c", "adt_CD14", "adt_CD45RA","adt_CD56","adt_CD8", "adt_CD16","adt_CD19","adt_CD34"))
plot_1 + theme(axis.text.x = element_text(angle = 45, hjust=1))
# We can compare the protein and mRNA levels of the same gene. Note that the name of genes could be different than the proteins they encode.
# Now we check the expression patterns of all 10 proteins in the clusters.
# In these plots, protein (ADT) levels are on top, and RNA levels are on the bottom
FeaturePlot(CBMC8K_human, reduction = "tsne",features = c("adt_CD3", "adt_CD11c", "adt_CD8", "adt_CD16", "HUMAN-CD3E", "HUMAN-ITGAX", "HUMAN-CD8A", 
                               "HUMAN-FCGR3A"), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)

FeaturePlot(CBMC8K_human, reduction = "tsne",features = c("adt_CD3", "adt_CD19", "adt_CD14", "adt_CD11c", "HUMAN-CD3E", "HUMAN-CD19", "HUMAN-CD14", 
                                                                     "HUMAN-ITGAX"), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)

FeaturePlot(CBMC8K_human, reduction = "tsne",features = c("adt_CD3", "adt_CD4", "adt_CD8", "HUMAN-CD3E", "HUMAN-CD4", "HUMAN-CD8A"
                                                                     ), min.cutoff = "q05", max.cutoff = "q95", ncol = 3)

FeaturePlot(CBMC8K_human, reduction = "tsne",features = c("adt_CD45RA", "adt_CD56", "adt_CD16", "adt_CD34", "HUMAN-PTPRC", "HUMAN-NCAM1", "HUMAN-FCGR3A", 
                                                                     "HUMAN-CD34"), min.cutoff = "q05", max.cutoff = "q95", ncol = 4)

# We can also use "RidgePlot" to visualize the expression profile of these proteins in each cell type.
RidgePlot(CBMC8K_human, features = c("adt_CD3", "HUMAN-CD3E"), ncol = 2)

# We can even draw scatter plots based on protein levels (like the classical biaxial plots for FACS). 
FeatureScatter(CBMC8K_human, feature1 = "adt_CD3", feature2 = "adt_CD19")
FeatureScatter(CBMC8K_human, feature1 = "adt_CD4", feature2 = "adt_CD8")

# view relationship between protein and RNA
FeatureScatter(CBMC8K_human, feature1 = "HUMAN-CD3E", feature2 = "adt_CD3")

# calculate the correlation
cor(CBMC8K_human[["RNA"]]@data["HUMAN-CD3E",], CBMC8K_human[["ADT"]]@data["CD3",])

# We can also restrict the pattern to one or more specific cell types:
tcells <- subset(CBMC8K_human, idents = c("Naive CD4 T", "Memory CD4 T", "CD8 T"))
FeatureScatter(tcells, feature1 = "HUMAN-CD3E", feature2 = "adt_CD3")
FeatureScatter(tcells, feature1 = "adt_CD4", feature2 = "adt_CD8")

# We can also draw a heatmap based on protein levels, similarly to what we did for the mRNA data.
CBMC8K_human_small <- subset(CBMC8K_human, downsample = 300)
adt.markers <- FindAllMarkers(CBMC8K_human_small, assay = "ADT", only.pos = TRUE)
DoHeatmap(CBMC8K_human_small, features = unique(adt.markers$gene), assay = "ADT", angle = 90) + NoLegend()

saveRDS(CBMC8K_human, file = "CBMC8K_human.rds")






# Cluster directly on protein levels #

# We clustered the single cells based on mRNA levels. We can also directly perform clustering on protein levels.

# To make it easier, we can switch the default assay to the 'ADT' assay, so that we don't need to specify it each time. 
DefaultAssay(CBMC8K_human) <- "ADT"

# We can run a PCA.
CBMC8K_human <- RunPCA(CBMC8K_human, features = rownames(CBMC8K_human), reduction.name = "pca_adt", reduction.key = "pca_adt_", 
               verbose = FALSE)
DimPlot(CBMC8K_human, reduction = "pca_adt")

# Because we only have 10 features here, instead of doing PCA, we can just use a standard euclidean distance matrix here.  
adt.data <- GetAssayData(CBMC8K_human, slot = "data")
adt.dist <- dist(t(adt.data))

# Before we re-cluster the data based on protein levels, we can store the currrent cluster IDs based on RNA expressions for later comparison.
CBMC8K_human[["rnaClusterID"]] <- Idents(CBMC8K_human)

# Using our distance matrix defined only on protein levels, we can now re-do tSNE using the following codes:
CBMC8K_human[["tsne_adt"]] <- RunTSNE(adt.dist, assay = "ADT", reduction.key = "adtTSNE_")
CBMC8K_human[["adt_snn"]] <- FindNeighbors(adt.dist)$snn
CBMC8K_human <- FindClusters(CBMC8K_human, resolution = 0.2, graph.name = "adt_snn")

# To annotate the protein clustering, we can compare the RNA and protein clustering:
clustering.table <- table(Idents(CBMC8K_human), CBMC8K_human$rnaClusterID)
clustering.table

# Based on the clustering table, we can now label the clusters:
new.cluster.ids <- c("CD4 T", "CD14+ Mono", "NK", "B", "CD8 T", "NK", "T/Mono doublets", "CD16+ Mono", "CD34+","pDCs", "B")
# Idents(CBMC8K_human) <- CBMC8K_human@meta.data$seurat_clusters
names(new.cluster.ids) <- levels(CBMC8K_human)
CBMC8K_human <- RenameIdents(CBMC8K_human, new.cluster.ids)
DimPlot(CBMC8K_human, reduction = "tsne_adt", label = TRUE, pt.size = 0.5) 

# Draw a heatmap using protein markers based on these clusters.
CBMC8K_human_small2 <- subset(CBMC8K_human, downsample = 300)
adt.markers2 <- FindAllMarkers(CBMC8K_human_small2, assay = "ADT", only.pos = TRUE)
DoHeatmap(CBMC8K_human_small2, features = unique(adt.markers2$gene), assay = "ADT", angle = 90) + NoLegend()

# Similar to what we did for the RNA clusters, we can visualize pretein levels based on the protein clusters:
tcells <- subset(CBMC8K_human, idents = c("CD4 T", "CD8 T"))
FeatureScatter(tcells, feature1 = "CD4", feature2 = "CD8")
RidgePlot(CBMC8K_human, features = c("adt_CD3", "adt_CD4", "adt_CD11c", "adt_CD14", "adt_CD45RA","adt_CD56","adt_CD8", "adt_CD16","adt_CD19","adt_CD34"), 
          ncol = 2)
RidgePlot(CBMC8K_human, features = c("adt_CD3", "adt_CD4"),ncol = 2)



# comparing RNA and protein clustering #

# At the end, let's compare the clustering based on either RNA or protein levels 
# of marker genes. Here we will visualize both the RNA and protein clustering based 
# on a tSNE generated using the ADT distance matrix.
tsne_rnaClusters <- DimPlot(CBMC8K_human, reduction = "tsne_adt", group.by = "rnaClusterID") +
  NoLegend()
tsne_rnaClusters <- tsne_rnaClusters + ggtitle("Clustering (scRNA-seq)") + 
  theme(plot.title = element_text(hjust = 0.5))
tsne_rnaClusters <- LabelClusters(plot = tsne_rnaClusters, id = "rnaClusterID", size = 4)

tsne_adtClusters <- DimPlot(CBMC8K_human, reduction = "tsne_adt", pt.size = 0.5) + 
  NoLegend()
tsne_adtClusters <- tsne_adtClusters + ggtitle("Clustering (ADT)") + 
  theme(plot.title = element_text(hjust = 0.5))
tsne_adtClusters <- LabelClusters(plot = tsne_adtClusters, id = "ident", size = 4)

CombinePlots(plots = list(tsne_rnaClusters, tsne_adtClusters), ncol = 2)

# Let reset the default assay back to "RNA" and save the Seurat object.
DefaultAssay(CBMC8K_human) <- "RNA"
saveRDS(CBMC8K_human, file = "CBMC8K_human.rds")

#                                                                     #
#                                                                     #
#######################################################################
