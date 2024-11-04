observeEvent(input$sample,{
  updateSelectInput(session,
                    "inputGene",
                    selected = c("Il4i1", "St7", "Isg15")
  )
})

observeEvent(input$clear,{
  updateSelectInput(session, "inputGene", selected = "")
})




# Submit Celltype
submit_btn_celltype <- observeEvent(input$submit_celltype,{
  allsigs <- read.xlsx("dataFiles/SuppTable3_CytokineSignatures.xlsx", input$cellType1)
  
  if (input$cellType1 == " ") {
    showNotification('Please select a cell type.', type = 'error')
  }
  else if (length(input$inputGene) == 0) {
    showNotification('Please enter in genes.', type = 'error')
  }
  else if (input$type_cytokines == 'inp' & length(input$inputCytokines) > 0){
    subset_genes <- subset(allsigs, Gene %in% input$inputGene)
    subset_data <- subset(subset_genes, Celltype_Str == input$cellType1)
    final_data <- subset(subset_data, Cytokine %in% input$inputCytokines)
    plot_height <- 500 + 2*nrow(final_data)
    
    p <- plot_ly(x=final_data$Cytokine, y=final_data$Gene, z = final_data$Avg_log2FC, type = "heatmap", height = plot_height) %>%
      layout(title = 'Heat Map Gene vs. Cytokine',
             yaxis = list(title = 'Gene'),
             xaxis = list(side = 'top',title = 'Cytokine'),
             margin = list(t=250, l=100))
    data$plot_input1 = p
    
    heatmap_data$selectedGenes = input$inputGene
  }
  else if (input$type_cytokines == 'inp' & length(input$inputCytokines) <= 0){
    showNotification('Please enter in cytokines.', type = 'error')
  }
  else if (input$type_cytokines == 'all'){
    subset_data <- subset(allsigs, Gene %in% input$inputGene)
    final_data <- subset(subset_data, Celltype_Str == input$cellType1)
    plot_height <- 500 + 2*nrow(final_data)
    p <- plot_ly(x=final_data$Cytokine, y=final_data$Gene, z = final_data$Avg_log2FC, type = "heatmap", height = plot_height) %>%
      layout(title = 'Heat Map Gene vs. Cytokine',
             yaxis = list(title = 'Gene'),
             xaxis = list(side = 'top',title = 'Cytokine'),
             margin = list(t=250, l=100))
    data$plot_all1 = p
    
    heatmap_data$selectedGenes = input$inputGene
  }
  else{
    return(NULL)
  }
})

# Submit Cytokines
submit_btn_cytokines <- observeEvent(input$submit_cytokines,{
  #allsigs <- read.xlsx("dataFiles/SuppTable3_CytokineSignatures.xlsx", "B_cell")
  allsigs <- dplyr::bind_rows(
    rio::import_list("dataFiles/SuppTable3_CytokineSignatures.xlsx", setclass = "tbl"),
    .id = "sheet"
  )
  
  if (input$cytokine1 == " ") {
    showNotification('Please select a cytokine.', type = 'error')
  }
  else if (length(input$inputGene) == 0) {
    showNotification('Please enter in genes.', type = 'error')
  }
  else if (input$type_cellType == 'inp' & length(input$inputCellTypes) > 0){
    subset_genes <- subset(allsigs, Gene %in% input$inputGene)
    subset_data <- subset(subset_genes, Cytokine == input$cytokine1)
    
    final_data <- subset(subset_data, Celltype_Str %in% input$inputCellTypes)
    plot_height <- 500 + 2*nrow(final_data)
    
    p <- plot_ly(x=final_data$Celltype_Str, y=final_data$Gene, z = final_data$Avg_log2FC, type = "heatmap", height = plot_height) %>%
      layout(title = 'Heat Map Gene vs. Cell Type',
             yaxis = list(title = 'Gene'),
             xaxis = list(side = 'top',title = 'Cell Type'),
             margin = list(t=250, l=100))
    data$plot_input2 = p
    
    heatmap_data$selectedGenes = input$inputGene
  }
  else if (input$type_cellType == 'all'){
    subset_data <- subset(allsigs, Gene %in% input$inputGene)
    final_data <- subset(subset_data, Cytokine == input$cytokine1)
    plot_height <- 500 + 2*nrow(final_data)
    p <- plot_ly(x=final_data$Celltype_Str, y=final_data$Gene, z = final_data$Avg_log2FC, type = "heatmap", height = plot_height) %>%
      layout(title = 'Heat Map Gene vs. Cell Type',
             yaxis = list(title = 'Gene'),
             xaxis = list(side = 'top',title = 'Cell Type'),
             margin = list(t=250, l=100))
    data$plot_all2 = p
    
    heatmap_data$selectedGenes = input$inputGene
  }
  else{
    return(NULL)
  }
})

# Taken out due to incompatibility with reactive value initialization
# Persists Selection Data
# observeEvent(input$type_cellType,{heatmap_data$selectedGenes = input$inputGene})
# observeEvent(input$type_cytokines,{heatmap_data$selectedGenes = input$inputGene})