observeEvent(input$sample,{
  updateTextAreaInput(session,
                      "inputGene",
                      value = "Ifitm3\nIsg15\nIfit3\nBst2\nSlfn5\nIsg20\nPhf11b\nZbp1\nRtp4"
  )
})

observeEvent(input$clear,{
  updateTextAreaInput(session, "inputGene", value = "")
})

# Submit Celltype
submit_btn_celltype <- observeEvent(input$submit_celltype,{
  if (input$cellType1 == " ") {
    showNotification('Please select a cell type.', type = 'error')
  }
  else if (input$inputGene == "") {
    showNotification('Please enter in genes.', type = 'error')
  }
  else if (input$type_cytokines == 'inp' & input$inputCytokines != ''){
    split_genes_vector <- split_input(input$inputGene)
    is_valid <- valid_comp(split_genes_vector, allsigs$gene)
    if (is_valid == 'invalid'){
      showNotification('Please enter in valid genes.', type = 'error')
    }
    else {
      if (is_valid != 'valid'){
        showNotification(HTML(paste0('The following genes were not found:', '<br>', paste(is_valid, collapse = ', '), '<br>', 'The calculation will proceed without these genes.')), type = 'error')
      }
      split_cytokines_vector <- split_input(input$inputCytokines)
      is_valid_cyto <- valid_comp(split_cytokines_vector, allsigs$sample)
      
      if (is_valid_cyto == 'invalid'){
        showNotification('Please enter in valid cytokines.', type = 'error')
      }
      else {
        if (is_valid_cyto != 'valid'){
          showNotification(HTML(paste0('The following cytokines were not found:', '<br>', paste(is_valid_cyto, collapse = ', '), '<br>', 'The calculation will proceed without these cytokines.')), type = 'error')
        }
        subset_genes <- subset(allsigs, gene %in% split_genes_vector)
        subset_data <- subset(subset_genes, celltype == input$cellType1)
        final_data <- subset(subset_data, sample %in% split_cytokines_vector)
        plot_height <- 500 + 2*nrow(final_data)
        
        p <- plot_ly(x=final_data$sample, y=final_data$gene, z = final_data$avg_logFC, type = "heatmap", height = plot_height) %>%
          layout(title = 'Heat Map Gene vs. Cytokine',
                 yaxis = list(title = 'Gene'),
                 xaxis = list(side = 'top',title = 'Cytokine'),
                 margin = list(t=250, l=100))
        data$plot_input1 = p
      }
    }
  }
  else if (input$type_cytokines == 'inp' & input$inputCytokines == ''){
    showNotification('Please enter in cytokines.', type = 'error')
  }
  else if (input$type_cytokines == 'all'){
    split_genes_vector <- split_input(input$inputGene)
    is_valid <- valid_comp(split_genes_vector, allsigs$gene)
    if (is_valid == 'invalid'){
      showNotification('Please enter in valid genes.', type = 'error')
    }
    else {
      if (is_valid != 'valid'){
        showNotification(HTML(paste0('The following genes were not found:', '<br>', paste(is_valid, collapse = ', '), '<br>', 'The calculation will proceed without these genes.')), type = 'error')
      }
      subset_data <- subset(allsigs, gene %in% split_genes_vector)
      final_data <- subset(subset_data, celltype == input$cellType1)
      plot_height <- 500 + 2*nrow(final_data)
      p <- plot_ly(x=final_data$sample, y=final_data$gene, z = final_data$avg_logFC, type = "heatmap", height = plot_height) %>%
        layout(title = 'Heat Map Gene vs. Cytokine',
               yaxis = list(title = 'Gene'),
               xaxis = list(side = 'top',title = 'Cytokine'),
               margin = list(t=250, l=100))
      data$plot_all1 = p
    }
    
  }
  else{
    return(NULL)
  }
})

# Submit Cytokines
submit_btn_cytokines <- observeEvent(input$submit_cytokines,{
  if (input$cytokine1 == " ") {
    showNotification('Please select a cytokine.', type = 'error')
  }
  else if (input$inputGene == "") {
    showNotification('Please enter in genes.', type = 'error')
  }
  else if (input$type_cellType == 'inp' & input$inputCellTypes != ""){
    split_genes_vector <- split_input(input$inputGene)
    is_valid <- valid_comp(split_genes_vector, allsigs$gene)
    if (is_valid == 'invalid'){
      showNotification('Please enter in valid genes.', type = 'error')
    }
    else {
      if (is_valid != 'valid'){
        showNotification(HTML(paste0('The following genes were not found:', '<br>', paste(is_valid, collapse = ', '), '<br>', 'The calculation will proceed without these genes.')), type = 'error')
      }
      split_cellTypes_vector <- split_input(input$inputCellTypes)
      is_valid_cell <- valid_comp(split_cellTypes_vector, allsigs$celltype)
      
      if (is_valid_cell == 'invalid'){
        showNotification('Please enter in valid cell types.', type = 'error')
      }
      else {
        if (is_valid_cell != 'valid'){
          showNotification(HTML(paste0('The following cell types were not found:', '<br>', paste(is_valid_cell, collapse = ', '), '<br>', 'The calculation will proceed without these cell types.')), type = 'error')
        }
        subset_genes <- subset(allsigs, gene %in% split_genes_vector)
        subset_data <- subset(subset_genes, sample == input$cytokine1)
        
        final_data <- subset(subset_data, celltype %in% split_cellTypes_vector)
        # print(final_data)
        plot_height <- 500 + 2*nrow(final_data)
        
        p <- plot_ly(x=final_data$celltype, y=final_data$gene, z = final_data$avg_logFC, type = "heatmap", height = plot_height) %>%
          layout(title = 'Heat Map Gene vs. Cell Type',
                 yaxis = list(title = 'Gene'),
                 xaxis = list(side = 'top',title = 'Cell Type'),
                 margin = list(t=250, l=100))
        #show('barid')
        data$plot_input2 = p
      }
    }
  }
  else if (input$type_cellType == 'all'){
    # split_cellTypes_vector <- c("B_cell", "cDC1", "cDC2", "eTAC", "ILC", "LEC", "Macrophage", "MigDC", "Monocyte", "Neutrophil", "NK_cell", "pDC", "T_cell_CD4", "T_cell_CD8", "T_cell_gd", "Treg")
    split_genes_vector <- split_input(input$inputGene)
    is_valid <- valid_comp(split_genes_vector, allsigs$gene)
    if (is_valid == 'invalid'){
      showNotification('Please enter in valid genes.', type = 'error')
    }
    else {
      if (is_valid != 'valid'){
        showNotification(HTML(paste0('The following genes were not found:', '<br>', paste(is_valid, collapse = ', '), '<br>', 'The calculation will proceed without these genes.')), type = 'error')
      }
      subset_data <- subset(allsigs, gene %in% split_genes_vector)
      final_data <- subset(subset_data, sample == input$cytokine1)
      plot_height <- 500 + 2*nrow(final_data)
      p <- plot_ly(x=final_data$celltype, y=final_data$gene, z = final_data$avg_logFC, type = "heatmap", height = plot_height) %>%
        layout(title = 'Heat Map Gene vs. Cell Type',
               yaxis = list(title = 'Gene'),
               xaxis = list(side = 'top',title = 'Cell Type'),
               margin = list(t=250, l=100))
      #show('barid')
      data$plot_all2 = p
    }
  }
  else{
    return(NULL)
  }
})