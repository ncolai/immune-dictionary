observeEvent(input$submit, {
  runjs("$('.busy').show();")
  hide("plot")
  hide("table")
  hide("download")
  hide("download_table")
}, priority = 1, ignoreInit = TRUE)

observeEvent(input$gene_file, {
  # if gene file is uploaded, display the data in the Text Area for submission
  table_data <- input$gene_file$datapath
  my_data <- read_excel(table_data)
  gene_column <- my_data %>% select(1)
  gene_vector <- gene_column[[1]]
  updateTextAreaInput(session, "inputGene", value = paste(gene_vector, collapse='\n'))
})

genes <- observeEvent(input$submit,{
  # hide tab B results
  hide("plot_B")
  hide("table_B")
  hide("download_B")
  hide("radio_btns_B")
  hide("download_table_B")
  
  if (input$inputGene != '' & input$inputCell != ' '){
      split_genes_vector <- split_input(input$inputGene)
      is_valid <- valid_genes(split_genes_vector)
      if (is_valid == 'invalid'){
        showNotification('Please enter in valid genes.', type = 'error')
      }
      else {
        if (is_valid != 'valid'){
          showNotification(HTML(paste0('The following genes were not found:', '<br>', paste(is_valid, collapse = ', '), '<br>', 'The calculation will proceed without these genes.')), type = 'error')
        }
        print(split_genes_vector)
        cat('Calculating Gene Set Enrichment Score...\n')
        tryCatch({
          if (input$inputMethod == 'Score'){
            # Use GeneSetEnrichmentScore for table and plot
            df_irea = GeneSetEnrichmentScore(split_genes_vector, input$inputCell)
            data$table = df_irea[,c(5,1,2,3)] %>%
              arrange(desc(ES)) %>%
              plyr::rename(c("ES" = "Enrichment Score"))
            data$compass_plot = IreaCompassPlot(df_irea)
          }
          else {
            # Use GeneSetEnrichmentHyperTest for table and plot
            df_irea = GeneSetEnrichmentHyperTest(split_genes_vector, input$inputCell)
            data$table = df_irea[,c(1,2,3,4)] %>%
              arrange(pval) %>%
              plyr::rename(c("ES" = "Enrichment Score"))
            data$compass_plot = IreaCompassPlot(df_irea, color_by = 'pval')
          }
        },
        error = function(e){
          showNotification('Error in calculation. Please try again with a different cell type.', duration = NULL, type = 'error')
          data$table = NULL
          data$compass_plot = NULL
        })
        
      }
  } 
  else {
    showNotification('Please enter in genes and select a cell type.', type = 'error')
  }
  # show end results
  show("plot")
  show("table")
  show("download")
  show("download_table")
  runjs("$('.busy').hide();")
}, ignoreInit = TRUE)

output$table <- renderDataTable({
  # display table data
  req(data$table)
  data$table
})

plotInput <- function(){
  # display plot data
  req(data$compass_plot)
  data$compass_plot
}

output$plot <- renderPlot({
    plotInput()
}, height = 500)

output$download <- renderUI({
  # if plot is shown, display download buttons for plot
  if(!is.null(plotInput())) {
    tagList(downloadButton('downloadPNG', 'Download as PNG'),
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

output$downloadPNG <- downloadHandler(
  filename = "irea_compass_plot.png",
  content = function(file) {
    png(file, width     = 6,
        height    = 6,
        units     = "in",
        res       = 1200,
        pointsize = 8)
    print(plotInput())
    dev.off()
  })    

output$downloadJPG <- downloadHandler(
  filename = "irea_compass_plot.jpg",
  content = function(file) {
    jpeg(file, width     = 6,
        height    = 6,
        units     = "in",
        res       = 1200,
        pointsize = 8)
    print(plotInput())
    dev.off()
  })

output$downloadTable <- downloadHandler(
  filename = function(){"df_irea.csv"}, 
  content = function(fname){
    write.csv(data$table, fname)
  }
)