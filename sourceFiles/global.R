# library(rdrop2)
# library(shinycssloaders)

# if (!exists('lig_seurat')){
#   cat('Loading data\n')
#   load('refdata.RData')
#   cat('Loaded data')
# }

if ((!exists('valid_gene_list')) || (!exists('valid_gene_list_human'))){
  load('dataFiles/valid_gene_list.RData')
}

lig_seurat_data <- vector(mode="list", length=15)
names(lig_seurat_data) <- readLines("sourceFiles/lig_seurat_data.txt")
