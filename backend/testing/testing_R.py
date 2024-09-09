import rpy2.robjects as robjects

r_script = """
library(Seurat)
library(SeuratDisk)
library(clustree)
set.seed(1234)

h5ad_path <- 'pbmc35.h5ad'
h5seurat_path <- 'pbmc35.h5seurat' # This is the corrected part
rds_path <- 'pbmc35.rds'

Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE) # Specify the output format and overwrite if exists

# Load the converted h5Seurat file as a Seurat object
seuratObject <- LoadH5Seurat(h5seurat_path) # Corrected to load the h5Seurat file



seuratObject <- NormalizeData(seuratObject)
seuratObject <- FindVariableFeatures(seuratObject)
seuratObject <- ScaleData(seuratObject)
seuratObject <- RunPCA(seuratObject)
seuratObject <- FindNeighbors(seuratObject)
seuratObject <- FindClusters(seuratObject, res = seq(0, 2, by = 0.1))
clustreePlot <- clustree(seuratObject, prefix = "RNA_snn_res.")

# Now, save the plot as an SVG file
ggsave("clustreePlot.svg", plot = clustreePlot, width = 10, height = 10)
"""

# Execute the R script
robjects.r(r_script)

print("Conversion from h5ad to RDS completed.")
