make_matrix <- function(df,rownames = NULL){
  # make matrix from data frame
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}

observeEvent(input$dropdown_btn,{
  # hide/show optional parameters on click
  toggle('genediff_cutoff')
})

observeEvent(input$submit_tabB, {
  # hide results during calculation
  runjs("$('.busy').show();")
  hide("plot_B")
  hide("table_B")
  hide("radio_btns_B")
  hide("download_B")
  hide("download_table_B")
}, priority = 1, ignoreInit = TRUE)

genesB <- observeEvent(input$submit_tabB, {
  `%notin%` = Negate(`%in%`)
  # hide tab A results
  hide("plot")
  hide("table")
  hide("download")
  hide("download_table")
  
  if ((!is.null(input$matrix_file) | input$sample_matrix) & input$inputCell_tabB != ' '){
    if (input$sample_matrix){
      # use sample matrix data file
      my_data <- read_excel('exampleFiles/gene_test.xlsx')
    }
    else{
      # use uploaded file
      table_data <- input$matrix_file$datapath
      my_data <- read_excel(table_data, .name_repair = ~ ifelse(nzchar(.x), .x, paste('Sample', seq_along(.x) - 1)))  # rename columns without names
    }
    cat('calculating')
    data$input_profile <- make_matrix(select(my_data,-1), pull(my_data,1))      # make matrix from input
    
    # check if file format is correct
    # make sure each column has a title and only numbers
    if (!isTRUE(all(is.numeric(data$input_profile) & -1 <= data$input_profile & data$input_profile <= 1))){
      data$input_profile <- NULL
      showNotification(HTML(paste0('Your data file format is incorrect.', '<br>', 'Please make sure the first column is genes, 
                                   and the following columns only contain numbers (in the range [-1 to 1]) pertaining to the sample.')), duration = NULL, type = 'error')
    }
    else{
    # check if genes are valid
    is_valid <- valid_genes(rownames(data$input_profile))
    
    if (is_valid == 'invalid'){
      showNotification('Please submit valid genes in the first column.', duration = NULL, type = 'error')
      data$input_profile <- NULL
      data$table_tabB = NULL
      data$compass_plot_B = NULL
    }
    else {
      # if there are invalid genes, show error message
      if (is_valid != 'valid'){
        showNotification(HTML(paste0('The following genes were not found:', '<br>', paste(is_valid, collapse = ', '), '<br>', 'The calculation will proceed without these genes.')), duration = 10, type = 'error')
      }
      
      tryCatch({
        # perform GetEnrichmentScoreProjection and IreaCompassPlot calculations
        data$table_tabB <- GetEnrichmentScoreProjection(data$input_profile, input$inputCell_tabB, genediff_cutoff = input$genediff_cutoff)
        df_irea_pd1 = subset(data$table_tabB, Sample == colnames(data$input_profile)[1])
        data$compass_plot_B <- IreaCompassPlot(df_irea_pd1)
      },
      error = function(e){
        showNotification('Error in calculation. Please try again with a different cell type.', duration = NULL, type = 'error')
        data$input_profile = NULL
        data$table_tabB = NULL
        data$compass_plot_B = NULL
      })
    }
    }
  }
  else{
    showNotification('Please upload gene matrix file and select a cell type.', type = 'error')
  }
  runjs("$('.busy').hide();")
  show("plot_B")
  show("table_B")
  show("radio_btns_B")
  show("download_B")
  show("download_table_B")
}, ignoreInit = TRUE)

output$table_B <- renderDataTable({
  # if table data has been calculated, display table
  req(data$table_tabB)
  data$table_tabB_subset <- subset(data$table_tabB, Sample == input$rb)
  data$table_tabB_subset
})

plotInput_B <- function(){
  # if table data has been calculated, display plot
  req(data$table_tabB)
  if (!is.null(input$rb)){
    df_irea_pd1 = subset(data$table_tabB, Sample == input$rb)
    data$compass_plot_B <- IreaCompassPlot(df_irea_pd1)
    data$compass_plot_B
  }
  else {
    data$compass_plot_B
  }
}

output$plot_B <- renderPlot({
  plotInput_B()
}, height = 500)

output$radio_btns_B <- renderUI({
  req(data$input_profile)
  radioButtons("rb", "Choose one:", choices = colnames(data$input_profile))   # default choice is first option
})

output$download_matrix <- downloadHandler(
  # download example gene_matrix file
  filename = function(){"ex_gene.xlsx"},
  content <- function(file) {
    file.copy("exampleFiles/gene_test.xlsx", file)
  },
  contentType = "application/xlsx"
)

# observeEvent(input$matrix_file, {
#   updateCheckboxInput(session, 'sample_matrix', value = 'FALSE')
# })

observeEvent(input$sample_matrix, {
  # session$sendCustomMessage("upload_txt", "SOME OTHER TEXT")
  if (input$sample_matrix == TRUE){
    hide('matrix_file')
  }
  else{
    show('matrix_file')
  }
})

# download button tab B
output$download_B <- renderUI({
  req(data$compass_plot_B)
  tagList(downloadButton('downloadPNG_B', 'Download as PNG'),
            downloadButton('downloadJPG_B', 'Download as JPG'),
          HTML("<br><br><br>"))
})

output$download_table_B <- renderUI({
  req(data$table_tabB, data$table_tabB_subset)
  tagList(
    downloadButton('downloadTable_B', "Download table"),
    downloadButton('downloadTable_B_all', "Download all samples"),
    HTML("<br><br><br>")
  )
})

output$downloadPNG_B <- downloadHandler(
  filename = "irea_compass_plot.png",
  content = function(file) {
    png(file, width     = 6,
        height    = 6,
        units     = "in",
        res       = 1200,
        pointsize = 8)
    print(plotInput_B())
    dev.off()
  })    

output$downloadJPG_B <- downloadHandler(
  filename = "irea_compass_plot.jpg",
  content = function(file) {
    jpeg(file, width     = 6,
         height    = 6,
         units     = "in",
         res       = 1200,
         pointsize = 8)
    print(plotInput_B())
    dev.off()
  })    

output$downloadTable_B <- downloadHandler(
  filename = function(){"df_irea.csv"},
  content = function(fname){
    write.csv(data$table_tabB_subset, fname)
  }
)

output$downloadTable_B_all <- downloadHandler(
  filename = function(){"df_irea_all.csv"},
  content = function(fname){
    write.csv(data$table_tabB, fname)
  }
)