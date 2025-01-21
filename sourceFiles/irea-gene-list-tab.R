observeEvent(input$submit_compass_list, {
  runjs("$('.busy').show();")
  hideElement("plot")
  hideElement("download")
  hideElement("table")
  hideElement("download_table")
}, priority = 1, ignoreInit = TRUE)

observeEvent(input$submit_radar_list, {
  runjs("$('.busy').show();")
  hideElement("plot")
  hideElement("download")
  hideElement("table")
  hideElement("download_table")
}, priority = 1, ignoreInit = TRUE)

observeEvent(input$dropdown_btn,{
  # hide/show optional parameters on click
  toggle('inputMethod')
})

observeEvent(input$gene_file, {
  # if gene file is uploaded, display the data in the Text Area for submission
  table_data <- input$gene_file$datapath
  my_data <- read_excel(table_data)
  gene_column <- my_data %>% select(1)
  gene_vector <- gene_column[[1]]
  updateTextAreaInput(session, "inputGene", value = paste(gene_vector, collapse='\n'))
})

genes <- observeEvent(input$submit_compass_list,{
  # hide tab B results
  hideElement("plot_B")
  hideElement("download_B")
  hideElement("table_B")
  hideElement("download_table_B")
  hideElement("radio_btns_B")
  
  #hide tab C results
  hideElement("plot_C")
  hideElement("download_C")
  hideElement("table_C")
  hideElement("download_table_C")
  
  if (input$inputGene != '' & input$inputCell != ' '){
    split_genes_vector <- split_input(input$inputGene)
    is_valid <- valid_genes(split_genes_vector, species = tolower(input$speciesInput))
    if (is_valid[1] == 'invalid'){
      showNotification(HTML(paste0('Please enter in valid genes.', '<br><br>', 'Example mouse genes: Isg15, Irf7, Ncr1', '<br>', 'Example human genes: ISG15, IRF7, NCR1')), type = 'error')
    }
    else {
      if (is_valid[1] != 'valid'){
        showNotification(HTML(paste0('The following genes were not found:', '<br>', paste(is_valid, collapse = ', '), '<br>', 'The calculation will proceed without these genes.', '<br><br>', 'Example mouse genes: Isg15, Irf7, Ncr1', '<br>', 'Example human genes: ISG15, IRF7, NCR1')), type = 'error')
      }
      print(split_genes_vector)
      cat('Calculating Gene Set Enrichment Score...\n')
      tryCatch({
        if (input$inputMethod == 'Score'){
          # Use GeneSetEnrichmentScore for table and plot
          df_irea = GeneSetEnrichmentScore(split_genes_vector, input$inputCell, species = tolower(input$speciesInput))
          data$table = df_irea[,c(5,1,2,3)] %>%
            arrange(desc(ES)) %>%
            plyr::rename(c("ES" = "Enrichment Score"))
          print(df_irea)
          data$irea_plot = IreaCompassPlot(df_irea, color_by = "pval")
        }
        else {
          # Use GeneSetEnrichmentHyperTest for table and plot
          df_irea = GeneSetEnrichmentHyperTest(split_genes_vector, input$inputCell, species = tolower(input$speciesInput))
          data$table = df_irea[,c(1,2,3,4)] %>%
            arrange(pval) %>%
            plyr::rename(c("ES" = "Enrichment Score"))
          data$irea_plot = IreaCompassPlot(df_irea, color_by = "pval")
        }
      },
      error = function(e){
        showNotification('Error in calculation. Please try again with a different cell type.', duration = NULL, type = 'error')
        data$table = NULL
        data$irea_plot = NULL
      })
      
    }
  } 
  else {
    showNotification('Please enter in genes and select a cell type.', type = 'error')
  }
  # show end results
  showElement("plot")
  showElement("download")
  showElement("table")
  showElement("download_table")
  runjs("$('.busy').hide();")
}, ignoreInit = TRUE)

genes <- observeEvent(input$submit_radar_list,{
  # hide tab B results
  hideElement("plot_B")
  hideElement("download_B")
  hideElement("table_B")
  hideElement("download_table_B")
  hideElement("radio_btns_B")
  
  #hide tab C results
  hideElement("plot_C")
  hideElement("download_C")
  hideElement("table_C")
  hideElement("download_table_C")
  
  if (input$inputGene != '' & input$inputCell != ' '){
    split_genes_vector <- split_input(input$inputGene)
    is_valid <- valid_genes(split_genes_vector, species = tolower(input$speciesInput))
    if (is_valid == 'invalid'){
      showNotification(HTML(paste0('Please enter in valid genes.', '<br><br>', 'Example mouse genes: Isg15, Irf7, Ncr1', '<br>', 'Example human genes: ISG15, IRF7, NCR1')), type = 'error')
    }
    else {
      if (is_valid != 'valid'){
        showNotification(HTML(paste0('The following genes were not found:', '<br>', paste(is_valid, collapse = ', '), '<br>', 'The calculation will proceed without these genes.', '<br><br>', 'Example mouse genes: Isg15, Irf7, Ncr1', '<br>', 'Example human genes: ISG15, IRF7, NCR1')), type = 'error')
      }
      print(split_genes_vector)
      cat('Calculating Gene Set Enrichment Score...\n')
      tryCatch({
        if (input$inputMethod == 'Score'){
          # Use GeneSetEnrichmentScore for table and plot
          df_irea = PolarizationGeneSetEnrichmentScore(split_genes_vector, input$inputCell, species = tolower(input$speciesInput))
          data$table = df_irea[,c(5,1,2,3)] %>%
            arrange(desc(ES)) %>%
            plyr::rename(c("ES" = "Enrichment Score"))
          print(df_irea)
          data$irea_plot = IreaRadarPlot(df_irea, input_celltype = input$inputCell)
        }
        else {
          # Use GeneSetEnrichmentHyperTest for table and plot
          df_irea = PolarizationGeneSetEnrichmentHyperTest(split_genes_vector, input$inputCell, species = tolower(input$speciesInput))
          data$table = df_irea[,c(1,2,3,4)] %>%
            arrange(pval) %>%
            plyr::rename(c("ES" = "Enrichment Score"))
          data$irea_plot = IreaRadarPlot(df_irea, input_celltype = input$inputCell)
        }
      },
      error = function(e){
        print(e)
        showNotification('Error in calculation. Please try again with a different cell type.', duration = NULL, type = 'error')
        data$table = NULL
        data$irea_plot = NULL
      })
      
    }
  } 
  else {
    showNotification('Please enter in genes and select a cell type.', type = 'error')
  }
  # show end results
  showElement("plot")
  showElement("download")
  showElement("table")
  showElement("download_table")
  runjs("$('.busy').hide();")
}, ignoreInit = TRUE)

output$table <- renderDataTable({
  # display table data
  req(data$table)
  data$table
})

plotInput <- function(){
  # display plot data
  req(data$irea_plot)
  data$irea_plot
}

output$plot <- renderPlot({
  plotInput()
}, height = 500)

output$download <- renderUI({
  # if plot is shown, display download buttons for plot
  if(!is.null(plotInput())) {
    tagList(downloadButton('downloadPDF', 'Download as PDF'),
            downloadButton('downloadPNG', 'Download as PNG'),
            downloadButton('downloadJPG', 'Download as JPG'),
            HTML("<br><br><br>"))
  }
})

output$download_table <- renderUI({
  # if table is shown, display table download button
  if(!is.null(data$table)) {
    tagList(
      downloadButton('downloadTable',"Download table"),
      HTML("<br><br><br>")
    )
  }
})


output$downloadPDF <- downloadHandler(
  filename = "IREA_plot.pdf",
  content = function(file) {
    cairo_pdf(file, 7, 6)
    print(plotInput())
    dev.off()
  })  


output$downloadPNG <- downloadHandler(
  filename = "IREA_plot.png",
  content = function(file) {
    png(file, width     = 7,
        height    = 6,
        units     = "in",
        res       = 1200,
        pointsize = 8)
    print(plotInput())
    dev.off()
  })    

output$downloadJPG <- downloadHandler(
  filename = "IREA_plot.jpg",
  content = function(file) {
    jpeg(file, width     = 7,
         height    = 6,
         units     = "in",
         res       = 1200,
         pointsize = 8)
    print(plotInput())
    dev.off()
  })

output$downloadTable <- downloadHandler(
  filename = function(){"IREA_output.csv"}, 
  content = function(fname){
    write.csv(data$table, fname)
  }
)
