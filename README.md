# PMID-32832598-lung-endothelial-cells
Extracting the endothelial cells from the PMID 32832598 lung single cell dataset.

This is done using the standard Seurat workflow (Stuart et al., 2019). To filter out endothelial cells, a list of endothelial markers was used (CLDN5, VWF, FLT1, PECAM1) as an input for SCINA (Zhang et al., 2019). The lung endothelial cell datasets is available as an rds file, with 4103 cells identified as being endothelial.
