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
      file_name = paste0('dataFiles/rdata-celltype/200417-ligands-alldata-seurat-p3-', input$umap_featureInput1, '.rds', sep = '')
      lig_seurat[[input_celltype]] <- readRDS(file_name)
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
    plot_df$color_use = ifelse(plot_df$sample == input$umap_featureInput2, "Selected", "Other")
    fig = ggplot(plot_df %>% arrange(color_use), aes(x = UMAP_1, y = UMAP_2)) + 
      geom_point(aes(color = color_use)) + 
      theme_nothing() + 
      theme(legend.position = "none") + 
      scale_color_manual(values = c("gray", "red")) +
      ggtitle(paste0(input$umap_featureInput1, " UMAP", sep = '')) +
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
      profiles_cc_gene = as.vector(lig_seurat[[input_celltype]]@assays[['RNA']][input$umap_gene_input,])
      
      plot_df = cbind(metadata_cc[, c("UMAP_1", "UMAP_2")], gene_exp = profiles_cc_gene)
      fig2 = ggplot(plot_df %>% arrange(gene_exp), aes(x = UMAP_1, y = UMAP_2)) + 
        geom_point(aes(color = gene_exp)) + 
        theme_nothing() + 
        theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold.italic")) + 
        scale_color_gradient(low = "gray", high = "blue") +
        ggtitle(input$umap_gene_input) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        labs(x = NULL, y = NULL)
      
      
      
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
      
      data$plot_res <- subplot(ggplotly(fig, height = 1200), 
                               ggplotly(fig2, height = 1200), 
                               nrows = 2, margin = 0.1) %>%
        layout(title = paste0("<b>", input$umap_featureInput1, " UMAP</b>", sep = ''), 
               annotations = annotations)
      
      # data$plot_res <- subplot(ggplotly(fig, height = 1200), 
      #                          ggplotly(fig2, height = 1200), 
      #                          nrows = 2, margin = 0.1) %>%
      #   layout(title = paste0("<b>", input$umap_featureInput1, " UMAP</b>", sep = ''), 
      #          annotations = annotations)
      umap_data$selectedGenes = input$umap_gene_input
    } else {
      data$plot_res <- ggplotly(fig, height = 600)
      umap_data$selectedGenes = input$umap_gene_input
    }
  }
})