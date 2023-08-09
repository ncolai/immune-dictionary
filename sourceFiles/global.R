# library(rdrop2)
# library(shinycssloaders)

# if (!exists('lig_seurat')){
#   cat('Loading data\n')
#   load('refdata.RData')
#   cat('Loaded data')
# }

if (!exists('valid_gene_list')){
  load('dataFiles/valid_gene_list.RData')
}

lig_seurat_data <- vector(mode="list", length=14)
names(lig_seurat_data) <- readLines("sourceFiles/lig_seurat_data.txt")

for (x in 1:14) {
  if (is.null(lig_seurat_data[[x]])){
    print('loading data')
    lig_seurat_data[[x]] <- readRDS(paste('dataFiles/ref_data/ref_data_', names(lig_seurat_data)[[x]], '.RDS', sep = ''))
  }
}

lig_seurat <- vector(mode="list", length=14)
names(lig_seurat) <- readLines("sourceFiles/lig_seurat_data.txt")

for (x in 1:14) {
  if (is.null(lig_seurat[[x]])){
    lig_seurat[[x]] <- readRDS(paste0("dataFiles/celltype/200417-ligands-alldata-seurat-p3-", 
                                                    names(lig_seurat)[[x]], ".rds", sep = ''))
  }
}