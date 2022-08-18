
#Set working directory
setwd("~/Documents/Shea Lab/Eiji's Spine Bridge Analysis /Illumina Files")

library(Seurat)
library(dplyr)
library(googledrive)
 #test test
# test 2

# get files from google drive 
#drive_download("")

# Input files from 10x
# reads output from Cell Ranger (pipelines process BCL files (seq output) into 
# FASTQ for 10x genomics method and can align reads to reference genome and 
# generate gene count matrix using STAR)
spleen.pbs.data <- Read10X(data.dir = "Spleen_PBS")
spleen.np.data <- Read10X(data.dir = "Spleen_NP")
bridge.pbs.data <- Read10X(data.dir = "Bridge_PBS")
bridge.np.data <- Read10X(data.dir = "Bridge_NP")

#Create Seurat objects
# use count matrix to create seurat object, which has data (ie count matrix)
# and analysis (ie PCA/clustering) 
# initialized with raw data 
spleen.pbs <- CreateSeuratObject(spleen.pbs.data, min.cells = 3)
spleen.pbs@meta.data$sample <- "spleen.pbs" # rename sample name
spleen.pbs <- NormalizeData(spleen.pbs) # normalize the data with log norm
spleen.pbs <- ScaleData(spleen.pbs, verbose = F) # linear transformation prior to
# dimensional reduction like PCA. Shifts expression so mean exp is 0 across all cells,
# scales exp, so var across cells is 1 
# verbose calls progress bar 

spleen.np <- CreateSeuratObject(spleen.np.data, min.cells = 3)
spleen.np@meta.data$sample <- "spleen.np"
spleen.np <- NormalizeData(spleen.np)
spleen.np <- ScaleData(spleen.np, verbose = F)

bridge.pbs <- CreateSeuratObject(bridge.pbs.data, min.cells = 3)
bridge.pbs@meta.data$sample <- "bridge.pbs"
bridge.pbs <- NormalizeData(bridge.pbs)
bridge.pbs <- ScaleData(bridge.pbs, verbose = F)

bridge.np <- CreateSeuratObject(bridge.np.data, min.cells = 3)
bridge.np@meta.data$sample <- "bridge.np"
bridge.np <- NormalizeData(bridge.np)
bridge.np <- ScaleData(bridge.np, verbose = F)

#Combine the two [FOUR?] Seurat objects and re-scale the data between the samples
# merges raw data count matrices of objects and makes new obj with combined raw
# count matrix. add.cell.ids appends given identifier so you can tell what original it is
data <- merge(x = spleen.pbs, y = c(spleen.np, bridge.pbs, bridge.np), 
              add.cell.ids = c('spleen.pbs', 'spleen.np', 'bridge.pbs', 'bridge.np'))

# Standard pre-processing workflow
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-") # percentage reads that map to MT genome
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # visualize QC metrics 
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10) # subset data
# filter out cells with unique features over 2500 or less than 200, filter cells >5% MT genes (probably dividing cell cycle)
data <- NormalizeData(data)
data <- ScaleData(data, verbose = F)

#Run PCA and tSNE on the data
data <- FindVariableFeatures(data) # calculate feature set with high cell to cell variation
data <- RunPCA(data, verbose = FALSE) # split into PCs based on variance (variable feat)
ElbowPlot(data) # rank PCs based on percentage of variance explained by each PC

# Cluster the cells
data <- FindNeighbors(data, dims = 1:9) 
data <- FindClusters(data, resolution = 0.5)

# Run non-linear dimension reduction, tSNE on same PCs as PCA 
data <- RunTSNE(data, check_duplicates = F) # group cells by similarity
TSNEPlot(data) # plot with all PCs
DotPlot(data, features = c("Ms4a1", "Itgax", "Siglecf", "Hbb-1", "Cd68", "Mrc1", "Ly6c2",
                           "Ly6g", "Ltf", "Retnlg", "Itga2", "Cd3e", "Cd4", "Cd8b1"))

#Save your Seurat object!!
save(data, file = 'np.sci.Robj')
# check if I can load it 
load("np.sci.Robj", verbose = TRUE)


# Finding differentially expressed features (cluster biomarkers)
#Identify the genes that define each cluster
all.markers <- FindAllMarkers(object = data, only.pos = TRUE, min.pct = 0.1, thresh.use = 0.25)
all.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.csv(all.markers, file = 'np.sci.all.markers.csv')
# find markers that define clusters via DE. 


#Assigning cell identities based on your research
data <- RenameIdents(data, '1' = 'B Cells', '9' = 'B Cells', '11' = 'Dendritic Cells', '13' = 
                       'Dendritic Cells (plasmacytoid)', '16' = 'Erythrocytes', '10' = 'Macrophages', 
                     '5' = 'Monocytes', '3' = 'Neutrophils', '8' = 'Neutrophils', '14' = 'Neutrophils',
                     '15' = 'Neutrophils', '7' = 'NK Cells', '4' = 'NK Cells (Ly6c+)', '17' = 'NK Cells (Ms4a1+)', 
                     '12' = 'NKT Cells', '0' = 'T Cells (Cd4+)', '6' = 'T Cells (Cd4+)', '2' = 'T Cells (Cd8+)')

TSNEPlot(data)

# Kate's Additions --------------------------------------------------------

# Feature Plots
#FeaturePlot(data), features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
#                               "CD8A"))

setwd("~/Documents/Shea Lab/Eiji's Spine Bridge Analysis ")
library(Seurat)
library(dplyr)

mac = subset(data, idents = "Macrophages")


load("mac.Robj", verbose = TRUE)
load("mon.Robj", verbose = TRUE)
load("neut.Robj", verbose = TRUE)

# Neutrophil feature plot
FeaturePlot(neut, features = c("Cd3e", "Ms4a1", "Ly6g", "C1qa", "Siglecf", "Fcer1a",
                               "Itga2", "Ltf", "Ifit3", "Cxcl2", "Mmp9", "Il1b"))
DotPlot(neut, features = c("Cd3e", "Ms4a1", "Ly6g", "C1qa", "Siglecf", "Fcer1a",
                           "Itga2", "Ltf", "Ifit3", "Cxcl2", "Mmp9", "Il1b"))
TSNEPlot(neut)
#UMAPPlot(neut)

# Macrophage feature plot
FeaturePlot(mac, features = c("Mrc1", "Cd86", "Arg1", "Tlr1", "Tlr8", "Vegfb", 
                              "Tlr2", "H2-Aa"))
DotPlot(mac, features = c("Mrc1", "Cd86", "Arg1", "Tlr1", "Tlr8", "Vegfb", 
                          "Tlr2", "H2-Aa", "Tnf", "Il1b", "Il6", "Nos2", "Ccr2", "Ccr7"))
prop.table(table(mac@active.ident, mac@meta.data$sample), margin = 2) # what two things you're comparing, margin = 2 normalizes across columns
# margin = 1 row normalized, no margin = total row & columns normalized 
# just table (without prop.table) will give regular counts 
barplot(prop.table(table(mac@active.ident, mac@meta.data$sample), margin = 2))
TSNEPlot(mac)

# Monocyte feature plot
FeaturePlot(mon, features = c("Ms4a1", "Cd3e", "Cx3cr1", "Ly6c2", "Cxcr2", "Ccr1"))
DotPlot(mon, features = c("Ms4a1", "Cd3e", "Cx3cr1", "Ly6c2", "Cxcr2", "Ccr1")) # easier to look at population markers 
TSNEPlot(mon)


# Neutrophil Analysis -----------------------------------------------------

library(UMAP)
# Neutrophil feature plot
FeaturePlot(neut, features = c("Cd3e", "Ms4a1", "Ly6g", "C1qa", "Siglecf", "Fcer1a",
                               "Itga2", "Ltf", "Ifit3", "Cxcl2", "Mmp9", "Il1b"))
DotPlot(neut, features = c("Cd3e", "Ms4a1", "Ly6g", "C1qa", "Siglecf", "Fcer1a",
                           "Itga2", "Ltf", "Ifit3", "Cxcl2", "Mmp9", "Il1b"))
TSNEPlot(neut)

NeutUmap <- RunUMAP(neut, dims = 1:10)
DimPlot(NeutUmap, reduction = "umap")
UMAPPlot(umap)


