library(ggplot2)
library(reshape2)
library(scales)
library(fst)
library(plotly)
library(qs)

cellList <- c(readLines("sourceFiles/lig_seurat_data.txt"))
cytokineList <- c(readLines("sourceFiles/cytokine_list.txt"), "")
geneList <- c(readLines("sourceFiles/gene_list.txt"))

# List of cell types to process
# cell_types <- c("B_cell", "cDC1", "cDC2", "eTAC", "Langerhans", "Macrophage", "Macrophage1", 
#                 "Macrophage2", "MigDC", "MigDC1", "Monocyte", "Neutrophil", "NK_cell", 
#                 "pDC", "T_cell_CD4", "T_cell_CD8", "T_cell_gd", "T_cell_gd1", "T_cell_gd2", "Treg")

cell_types <- c("B_cell", "cDC1", "cDC2", "Macrophage",
                "MigDC", "Monocyte", "Neutrophil", "NK_cell",
                "pDC", "T_cell_CD4", "T_cell_CD8", "T_cell_gd", "Treg")

cell_types <- c("MigDC")

# cell_types <- c("B_cell")
# violinInputGene <- c("Il4i1", "St7", "Isg15")
base_data_dir <- "dataFiles/violin_plot_creation/"
dir.create(base_data_dir, recursive = TRUE, showWarnings = FALSE)

for (cell_type in cell_types) {
  cell_type_dir <- file.path(base_data_dir, cell_type)
  dir.create(cell_type_dir, recursive = TRUE, showWarnings = FALSE) 
  
  fst_filename <- paste0("dataFiles/celltype_fst/200417-ligands-alldata-seurat-p3-", cell_type, ".rds.fst")
  if (file.exists(fst_filename)) {
    cell_data <- read_fst(fst_filename)
    for (input_violinInputGene in geneList){
    plot_df <- cell_data[c(input_violinInputGene, 'celltype', 'sample', 'rep')]
    plot_df_melt <- reshape2::melt(plot_df, id.vars = names(subset(plot_df, select = c('celltype','sample','rep'))))
     
    violin_plot_results <- lapply(input_violinInputGene, function(myGene) {
      plot_filename <- file.path(cell_type_dir, paste0(cell_type, "_violin_plot_", myGene, ".qs"))
      
      if (!file.exists(plot_filename)) {
        filtered_df_all <- subset(plot_df_melt, variable == myGene)
        # filtered_df_all[filtered_df_all<0] = 0
        color_palette <- hue_pal()(length(cytokineList))
        p <- ggplot(filtered_df_all, aes(x = sample, y = value, fill = sample)) +
          geom_violin(scale = "width", linewidth = 0) +
          theme_minimal() +
          theme(axis.text.x = element_text(size = 4.5, angle = 90, vjust = 0.5, hjust=1),
                axis.text.y = element_text(size = 6),
                axis.title.x = element_blank(),
                axis.title.y = element_text(size = 7, color = "black"),
                strip.text = element_text(face = "italic", color = "black", size = 5.5),
                strip.placement = "inside",
                strip.background = element_blank(),
                panel.spacing = unit(2, "mm"),
                legend.position = "none") +
          scale_fill_manual(values = color_palette) +
          facet_grid(~variable, scale = "free_y") +
          scale_y_continuous(breaks = c(0, 3, 6)) +
          ylab("Expression")
        fig <- ggplotly(p)
        
        qs::qsave(fig, plot_filename)
        # print(fig)
      }
    })
    }
  }
}


#testing
# (violinplot_test <- qread(file.path("dataFiles/violin_plot_creation/B_cell", paste0("B_cell", "_violin_plot_", "St7", ".qs"))))
