observeEvent(input$violinSample,{
  updateTextAreaInput(session,
                      "violinInputGene",
                      value = "Ifitm3\nIsg15\nIfit3\nBst2\nSlfn5\nIsg20\nPhf11b\nZbp1\nRtp4"
  )
})

observeEvent(input$violinClear,{
  updateTextAreaInput(session, "violinInputGene", value = "")
})

output$download_violin_all <- renderUI({
  req(data$violin_plot_all)
  tagList(downloadButton('download_violins_all_csv', 'Download Violin Plot as CSV'),
          downloadButton('download_violins_all_pdf', 'Download Violin Plot as PDF'),
          HTML("<br><br><br>"))
})

output$download_violin_inp <- renderUI({
  req(data$violin_plot_input)
  tagList(downloadButton('download_violins_inp_csv', 'Download Violin Plot as CSV'),
          downloadButton('download_violins_inp_pdf', 'Download Violin Plot as PDF'),
          HTML("<br><br><br>"))
})

output$download_violins_all_csv <- downloadHandler(
  filename = function(){"violin_plot_all.csv"},
  content = function(fname){
    write.csv(data$violin_all_data, fname)
  }
)

output$download_violins_inp_csv <- downloadHandler(
  filename = function(){"violin_plot_inp.csv"},
  content = function(fname){
    write.csv(data$violin_input_data, fname)
  }
)

final_plot_violin <- function(num) {
  if (num == 1) {
    req(data$violin_plot_all)
    plot_all <- data$violin_plot_all
    plot_all
  } else if (num == 2) {
    req(data$violin_plot_inp)
    plot_inp <- data$violin_plot_inp
    plot_inp
  }
}

output$download_violins_all_pdf <- downloadHandler(
  filename = "violin_plot_all.pdf",
  content = function(file) {
    plot_to_save <- final_plot_violin(1)
    save_image(plot_to_save, file = file, scale = 5)
  }) 

output$download_violins_inp_pdf <- downloadHandler(
  filename = "violin_plot_inp.pdf",
  content = function(file) {
    plot_to_save <- final_plot_violin(2)
    save_image(plot_to_save, file = file, scale = 5)
  }) 

submit_btn_celltype2 <- observeEvent(input$violinSubmit,{
  if (input$cellType2 == " ") {
    showNotification('Please select a cell type.', type = 'error')
  }
  else if (input$violinInputGene == "") {
    showNotification('Please enter in genes.', type = 'error')
  }
  else if (input$violin_type_cytokines == 'inp' & input$violin_inputCytokines != ''){
    
    if (is.null(lig_seurat[[input$cellType2]])){
      lig_seurat[[input$cellType2]] <- readRDS(paste0("dataFiles/celltype/200417-ligands-alldata-seurat-p3-", 
                                                      input$cellType2, ".rds", sep = ''))
    }
    
    split_genes_vector <- split_input(input$violinInputGene)
    is_valid <- valid_comp(split_genes_vector, allsigs$gene)
    if (is_valid == 'invalid'){
      showNotification('Please enter in valid genes.', type = 'error')
    }
    else {
      if (is_valid != 'valid'){
        showNotification(HTML(paste0('The following genes were not found:', '<br>', paste(is_valid, collapse = ', '), '<br>', 'The calculation will proceed without these genes.')), type = 'error')
      }
      split_cytokines_vector <- split_input(input$violin_inputCytokines)
      is_valid_cyto <- valid_comp(split_cytokines_vector, allsigs$sample)
      
      if (is_valid_cyto == 'invalid'){
        showNotification('Please enter in valid cytokines.', type = 'error')
      }
      else {
        if (is_valid_cyto != 'valid'){
          showNotification(HTML(paste0('The following cytokines were not found:', '<br>', paste(is_valid_cyto, collapse = ', '), '<br>', 'The calculation will proceed without these cytokines.')), type = 'error')
        }
        
        raw_data = lig_seurat[[input$cellType2]]
        plot_df = cbind(t(as.matrix(raw_data@assays[['RNA']][split_genes_vector,])),
                        raw_data@meta.data)
        plot_df_melt = melt(plot_df, id.vars = names(raw_data@meta.data))
        
        data$violin_input_data <- plot_df_melt[c('celltype', 'sample', 'rep', 'variable', 'value')]
        
        final_data <- list()
        
        for (i in 1:length(split_genes_vector)) {
          sub <- subset(plot_df_melt, variable == split_genes_vector[[i]])
          filtered_df_inp <- sub[sub$sample %in% split_cytokines_vector,]
          
          fig <- 
            plot_ly(x=filtered_df_inp$sample, 
                    y=filtered_df_inp$value,
                    split=filtered_df_inp$sample,
                    type = "violin",
                    box = list(visible = T),
                    meanline = list(visible = T),
                    height = 100*length(split_genes_vector),
                    marker = list(symbol = 'line-ns'),
                    # scalemode = 'width',
                    line = list(width = 0)) %>%
            layout(title = paste0('Violin Plot Gene Expression vs. Cytokine'),
                   xaxis = list(side = 'bottom', title = '', tickangle = 45),
                   yaxis = list(rangemode = 'nonnegative'),
                   showlegend = FALSE,
                   font=list(size=6),
                   annotations = list(x = 0.5,
                                      y = 0.9,
                                      xref = "paper",
                                      yref = "paper",
                                      xanchor = "center",
                                      yanchor = "bottom",
                                      text = split_genes_vector[[i]],
                                      showarrow = FALSE))
          
          fig2 <- plot_ly()%>%
            layout(xaxis=list(visible="FALSE",
                              color="rgba(0,0,0,0)",
                              tickfont =list(color="rgba(0,0,0,0)"),
                              showgrid = F),
                   yaxis=list(visible="FALSE",
                              color="rgba(0,0,0,0)",
                              tickfont =list(color="rgba(0,0,0,0)"),
                              showgrid = F),
                   showlegend = FALSE)
          
          fig3 <- subplot(fig, fig2, nrows = 2, heights = c(0.8, 0.2))
          
          final_data[[i]] <- fig3
        }
      
      data$violin_plot_input = subplot(final_data, 
                                       nrows = length(split_genes_vector))
      
      showNotification(HTML(paste0(
        'To zoom, select the section you would like to view in detail.')), 
        type = 'message')
      }
    }
  }
  else if (input$violin_type_cytokines == 'inp' & input$violin_inputCytokines == ''){
    showNotification('Please enter in cytokines.', type = 'error')
  }
  else if (input$violin_type_cytokines == 'all'){
    
    if (is.null(lig_seurat[[input$cellType2]])){
      lig_seurat[[input$cellType2]] <- readRDS(paste0("dataFiles/celltype/200417-ligands-alldata-seurat-p3-", 
                                                      input$cellType2, ".rds", sep = ''))
    }
    
    split_genes_vector <- split_input(input$violinInputGene)
    is_valid <- valid_comp(split_genes_vector, allsigs$gene)
    if (is_valid == 'invalid'){
      showNotification('Please enter in valid genes.', type = 'error')
    }
    else {
      if (is_valid != 'valid'){
        showNotification(HTML(paste0('The following genes were not found:', '<br>', paste(is_valid, collapse = ', '), '<br>', 'The calculation will proceed without these genes.')), type = 'error')
      }
      
      raw_data = lig_seurat[[input$cellType2]]
      plot_df = cbind(t(as.matrix(raw_data@assays[['RNA']][split_genes_vector,])),
                      raw_data@meta.data)
      plot_df_melt = melt(plot_df, id.vars = names(raw_data@meta.data))
      
      data$violin_all_data <- plot_df_melt[c('celltype', 'sample', 'rep', 'variable', 'value')]
      
      final_data <- list()
      
      for (i in 1:length(split_genes_vector)) {
        filtered_df_all <- subset(plot_df_melt, variable == split_genes_vector[[i]])
        
        fig <- 
          plot_ly(x=filtered_df_all$sample, 
                y=filtered_df_all$value,
                split=filtered_df_all$sample,
                type = "violin",
                box = list(visible = T),
                meanline = list(visible = T),
                height = 200*length(split_genes_vector),
                marker = list(symbol = 'line-ns'),
                scalemode = 'count',
                line = list(width = 0),
                spanmode = 'hard') %>%
          layout(title = paste0('Violin Plot Gene Expression vs. Cytokine'),
                xaxis = list(side = 'bottom', title = '', tickangle = 45),
                yaxis = list(rangemode = 'nonnegative'),
                showlegend = FALSE,
                font=list(size=6),
                annotations = list(x = 0.5,
                                   y = 0.9,
                                   xref = "paper",
                                   yref = "paper",
                                   xanchor = "center",
                                   yanchor = "bottom",
                                   text = split_genes_vector[[i]],
                                   showarrow = FALSE))
        
        fig2 <- plot_ly()%>%
          layout(xaxis=list(visible="FALSE",
                            color="rgba(0,0,0,0)",
                            tickfont =list(color="rgba(0,0,0,0)"),
                            showgrid = F),
                 yaxis=list(visible="FALSE",
                            color="rgba(0,0,0,0)",
                            tickfont =list(color="rgba(0,0,0,0)"),
                            showgrid = F),
                 showlegend = FALSE)

        fig3 <- subplot(fig, fig2, nrows = 2, heights = c(0.8, 0.2))
        
        final_data[[i]] <- fig3
      }
      
      data$violin_plot_all = subplot(final_data,
                                     nrows = length(split_genes_vector))

      showNotification(HTML(paste0(
        'To zoom, select the section you would like to view in detail.')), 
        type = 'message')
    }
    
  }
  else{
    return(NULL)
  }
})