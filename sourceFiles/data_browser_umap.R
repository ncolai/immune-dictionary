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
    start_time <- Sys.time()
    input_celltype = input$umap_featureInput1
    #logjs(paste('This is a test message!!!!', input_celltype))
    if (is.null(lig_seurat[[input_celltype]])){
      #print('loading data')
      #file_name = paste0('dataFiles/celltype/200417-ligands-alldata-seurat-p3-', input$umap_featureInput1, '.rds', sep = '')
      #lig_seurat[[input_celltype]] <- readRDS(file_name)
      file_name <- paste0("dataFiles/celltypeSmallUmap/",input$umap_featureInput1,".qs")
      lig_seurat[[input_celltype]] <- qread(file_name)
      Idents(lig_seurat[[input_celltype]]) <- "cluster" #fix docker img issue
      #print('data is loaded!!!')
    }
    
    if (input$umap_featureInput2 == 'ALL'){
      baseplot <- DimPlot(lig_seurat[[input_celltype]], reduction = "umap", pt.size = 1, group.by = "sample", label = FALSE)
    }
    else{
      names <- readLines("sourceFiles/cytokine_list.txt")
      values <- rep('red', times = length(names))
      inp_cytokine <- list()
      
      for (i in seq_along(names)) {
        inp_cytokine[[names[[i]]]] <- setNames(values[i], names[i])
      }
      
      baseplot <- DimPlot(lig_seurat[[input_celltype]], reduction = "umap", pt.size = 1, group.by = "sample", cols = inp_cytokine[[input$umap_featureInput2]], order = input$umap_featureInput2, label = FALSE)
    }
    fig <- baseplot + ggtitle(paste0(input$umap_featureInput1, " UMAP", sep = '')) + NoLegend() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                                                                        axis.title.x=element_blank(),
                                                                                                        axis.title.y=element_blank(),legend.position="none",
                                                                                                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                                                                                        panel.grid.minor=element_blank(),plot.background=element_blank())
    
    if (input$umap_gene_input != '' & valid_comp(input$umap_gene_input, geneList)){
      # check if gene is valid
      fig2 <- FeaturePlot(lig_seurat[[input_celltype]], features = c(input$umap_gene_input), order = TRUE, pt.size = 1) + NoLegend() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                                                             axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                                                                                                             axis.title.x=element_blank(),
                                                                                                                                             axis.title.y=element_blank(),legend.position="none",
                                                                                                                                             panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                                                                                                                             panel.grid.minor=element_blank(),plot.background=element_blank(),
                                                                                                                                             plot.margin = unit(c(0,0,0,0), "cm"))
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
      
      data$plot_res <- subplot(ggplotly(fig + theme(plot.margin = unit(c(0,3,0,3), "cm")), height = 900), 
                               ggplotly(fig2 + theme(plot.margin = unit(c(0,3,0,3), "cm")), height = 900), 
                               nrows = 2, margin = 0.05) %>%
        layout(title = paste0("<b>", input$umap_featureInput1, " UMAP</b>", sep = ''), 
               annotations = annotations)
      umap_data$selectedGenes = input$umap_gene_input
    } else {
      data$plot_res <- ggplotly(fig, height = 600)
      umap_data$selectedGenes = input$umap_gene_input
    }
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    print(elapsed_time)
    }
})
