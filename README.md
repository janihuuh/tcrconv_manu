# Determining recognition between TCRs and epitopes with transformer models 

Code related to scRNA+TCRab-seq analysis in TCRconv manuscript with the Liao et al data (GEO GSE145926). 

The scRNA-seq data is analyzed mainly with Python package scVI tools (v 0.14.5) (Gayoso et al., 2022) and R package Seurat (v 4.0.4) (Stuart et al., 2019). 

Cells with > 10 % mitochondrial gene counts, < 1,000 UMI counts, < 200 or > 6,000 detected genes, and cells with no detected TCR were filtered out. The highly variable genes were identified with "highly_variable_genes" function in scVI tools with default parameters, which were then used to learn latent embeddings with "model.SCVI" function in scVI tools with default parameters. 

The CD8+ T cells were then identified with SingleR (v 1.6.1) (Aran et al., 2019) and clustering, and the process was repeated with scVI tools. The gained embeddings were then used for finding clusters with "FindNeighbors" and "FindClusters" functions and further visualized with UMAP dimensionality reduction with "RunUMAP" function using default parameters in Seurat. 

The final optimal clustering threshold was chosen as 0.2 based on visual inspection of the clustering results in the UMAP reduced space. The markers used to define the clusters were found with Student’s t-test using the "FindMarkers" function in Seurat with logfc.threshold = 0.25 from experssion data that was scaled with "ScaleData" function with scaling factor of 10,000. 
