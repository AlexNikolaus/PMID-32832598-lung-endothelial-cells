library(Seurat)
library(tidyverse)
library(Matrix)
library(SCINA)

#Extract controls from the metadata file
table(GSE135893_IPF_metadata$Diagnosis)
Patients <- split(GSE135893_IPF_metadata, GSE135893_IPF_metadata$Diagnosis)
list2env(Patients, envir = .GlobalEnv)
table(Control$orig.ident)
table(GSE135893_IPF_metadata$orig.ident)
Controls <- unique(Control$orig.ident)

#Create object & subset for control
Lungs <- Read10X(data.dir = 'D:/Endothelial scRNA-Seq project/PMID 32832598 lungs', gene.column = 1)
Lungs <- CreateSeuratObject(Lungs, min.cells = 3, min.features = 200, project = 'Lungs')
Lungs <- subset(Lungs, idents = Controls)

#Standard Seurat workflow
Lungs <- NormalizeData(Lungs, normalization.method = "LogNormalize", scale.factor = 10000)
Lungs <- FindVariableFeatures(Lungs, selection.method = "vst", nfeatures = 2000)
Lungs <- ScaleData(Lungs)
Lungs <- RunPCA(Lungs, features = VariableFeatures(object = Lungs))
Lungs <- FindNeighbors(Lungs, dims = 1:10)
Lungs <- FindClusters(Lungs, resolution = 0.1)
Lungs <- RunUMAP(Lungs, dims = 1:10)

#Subset the clusters which have high marker gene expression
DimPlot(Lungs, reduction = "umap", label = TRUE)
VlnPlot(Lungs, features = c('CLDN5', 'VWF', 'FLT1', 'PECAM1'), ncol = 2)
Lungs <- subset(Lungs, idents = '3')

#Create list of markers for validation
Endothelial_Markers <- list(c('CLDN5', 'VWF', 'FLT1', 'PECAM1'))
names(Endothelial_Markers) <- 'Endothelial cells'

#Run SCINA to predict cell type
SCINA_results <- SCINA(Lungs@assays$RNA@data,
                       Endothelial_Markers,
                       max_iter = 2000, 
                       convergence_n = 100, 
                       convergence_rate = 0.999, 
                       sensitivity_cutoff = 0.9, 
                       rm_overlap=FALSE, 
                       allow_unknown=TRUE)
Lungs$cell_labels <- SCINA_results$cell_labels
DimPlot(Lungs,reduction = "umap", pt.size = 1, label = TRUE, group.by = 'cell_labels')

#Subset for endothelial cells
Lungs <- subset(Lungs, cell_labels == 'Endothelial cells')
Lungs$organ <- 'Lungs'
write_rds(Lungs, 'Lungs.rds')
