valid_comp <- function(gene, valid_genes) {
  if (gene %in% valid_genes){
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

observeEvent(input$umap_dropdown_btn,{
  # hide/show optional parameters on click
  toggle('umap_gene_input')
})

observeEvent(input$umap_submit, {
  if (input$umap_featureInput1 != ' ' & input$umap_featureInput2 != ' '){
    input_celltype = input$umap_featureInput1
    if (is.null(lig_seurat[[input_celltype]])){
      print('loading data')
      file_name = paste0("dataFiles/celltypeSmallUmap/",input$umap_featureInput1,".qs")
      # file_name = paste0('dataFiles/rdata-celltype/200417-ligands-seurat-', input$umap_featureInput1, '.RData', sep = '')
      lig_seurat[[input_celltype]] <- qread(file_name)
      print('data loaded')
      metadata_cc = cbind(lig_seurat[[input_celltype]]@meta.data, lig_seurat[[input_celltype]]@reductions$umap@cell.embeddings)
      
    }
    
    # if (input$umap_featureInput2 == 'ALL'){
    #   baseplot <- DimPlot(lig_seurat[[input_celltype]], reduction = "umap", pt.size = 1, group.by = "sample", label = FALSE)
    # }
    # else{
    # names <- readLines("sourceFiles/cytokine_list.txt")
    # values <- rep('red', times = length(names))
    # inp_cytokine <- list()
    # 
    # for (i in seq_along(names)) {
    #   inp_cytokine[[names[[i]]]] <- setNames(values[i], names[i])
    # }
    # 
    # baseplot <- DimPlot(lig_seurat[[input_celltype]], reduction = "umap", pt.size = 1, group.by = "sample", cols = inp_cytokine[[input$umap_featureInput2]], order = input$umap_featureInput2, label = FALSE)
    # }
    # fig <- baseplot + ggtitle(paste0(input$umap_featureInput1, " UMAP", sep = '')) + NoLegend() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
    #                                                                                                             axis.text.y=element_blank(),axis.ticks=element_blank(),
    #                                                                                                             axis.title.x=element_blank(),
    #                                                                                                             axis.title.y=element_blank(),legend.position="none",
    #                                                                                                             panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
    #                                                                                                             panel.grid.minor=element_blank(),plot.background=element_blank())
    
    plot_df = metadata_cc
    plot_df$color_use = ifelse(plot_df$sample == input$umap_featureInput2, input$umap_featureInput2, plot_df$sample)
    fig = ggplot(plot_df %>% arrange(sample), aes(x = UMAP_1, y = UMAP_2)) + 
      geom_point(aes(color = ifelse(sample == input$umap_featureInput2, "red", "#7F7F7F"),
                     text = paste0("UMAP_1: ", UMAP_1, 
                                   "\nUMAP_2: ", UMAP_2, 
                                   "\nSample: ", ifelse(sample == input$umap_featureInput2, input$umap_featureInput2, sample))),
                 size = 0.7) + 
      theme_nothing() + 
      theme(legend.position = "none",
            plot.margin = margin(t = 40, r = 0, b = 0, l = 0)) + 
      scale_color_manual(values = c("#7F7F7F", "red")) +
      ggtitle(paste0("<b>", input$umap_featureInput1, " UMAP</b>", sep = '')) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      labs(x = NULL, y = NULL)
    #fig  
    
    
    if (input$umap_gene_input != '' & valid_comp(input$umap_gene_input, geneList)){
      # check if gene is valid
      # fig2 <- FeaturePlot(lig_seurat[[input_celltype]], features = c(input$umap_gene_input)) + NoLegend() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
      #                                                                                                             axis.text.y=element_blank(),axis.ticks=element_blank(),
      #                                                                                                             axis.title.x=element_blank(),
      #                                                                                                             axis.title.y=element_blank(),legend.position="none",
      #                                                                                                             panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
      #                                                                                                             panel.grid.minor=element_blank(),plot.background=element_blank())
      # profiles_cc_gene = as.vector(lig_seurat[[input_celltype]]@assays[['RNA']][input$umap_gene_input,])
      # 
      # plot_df = cbind(metadata_cc[, c("UMAP_1", "UMAP_2")], gene_exp = profiles_cc_gene)
      # fig2 = ggplot(plot_df %>% arrange(gene_exp), aes(x = UMAP_1, y = UMAP_2)) + 
      #   geom_point(aes(color = gene_exp)) + 
      #   theme_nothing() + 
      #   theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold.italic")) + 
      #   scale_color_gradient(low = "gray", high = "blue") +
      #   ggtitle(input$umap_gene_input) +
      #   scale_x_continuous(expand=c(0,0)) +
      #   scale_y_continuous(expand=c(0,0)) +
      #   labs(x = NULL, y = NULL)
      
      umap_coordinates <- lig_seurat[[input_celltype]]@reductions$umap@cell.embeddings
      
      # Extract gene expression data
      expression_matrix <- lig_seurat[[input_celltype]]@assays[["RNA"]]@data
      if (input$umap_gene_input %in% rownames(expression_matrix)) {
        gene_expression <- expression_matrix[input$umap_gene_input, ]
      } else {
        stop(paste("Gene", input$umap_gene_input, "not found in the expression matrix"))
      }
      
      # Combine UMAP coordinates and gene expression into a data frame
      plot_data <- data.frame(UMAP1 = umap_coordinates[, 1],
                              UMAP2 = umap_coordinates[, 2],
                              Expression = gene_expression)
      
      plot_data <- plot_data[order(plot_data$Expression), ]
      
      if (all(plot_data$Expression == 0)) {
        color_scale <- scale_color_gradient(low = "#D3D3D3", high = "#D3D3D3")
      } else {
        color_scale <- scale_color_gradient(low = "#D3D3D3", high = "purple")
      }
      
      # Create the plot using ggplot2
      fig2 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = Expression)) +
        geom_point(size = 1) +
        color_scale +
        theme_minimal() +
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none",
              panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_blank())
      # plot.margin = unit(c(0, 0, 0, 0), "cm"))
      
      
      
      annotations = list( 
        list( 
          x = 0.5,  
          y = 0.5,  
          text = paste0("<b>", input$umap_gene_input, "</b>", sep = ''),  
          xref = "paper",  
          yref = "paper",  
          xanchor = "center",  
          yanchor = "bottom",  
          showarrow = FALSE 
        )
      )
      
      data$plot_res <- subplot(ggplotly(fig, tooltip = "text"), 
                               ggplotly(fig2), 
                               nrows = 2) %>%
        layout(title = paste0("<b>", input$umap_featureInput1, " UMAP</b>", sep = ''), 
               annotations = annotations,
               margin = list(t = 100),
               height = 1050)
      
      # data$plot_res <- subplot(ggplotly(fig, height = 1200), 
      #                          ggplotly(fig2, height = 1200), 
      #                          nrows = 2, margin = 0.1) %>%
      #   layout(title = paste0("<b>", input$umap_featureInput1, " UMAP</b>", sep = ''), 
      #          annotations = annotations)
      umap_data$selectedGenes = input$umap_gene_input
    } else {
      data$plot_res <- ggplotly(fig, tooltip = "text") %>%
        layout(height = 575)
      umap_data$selectedGenes = input$umap_gene_input
    }
  }
})