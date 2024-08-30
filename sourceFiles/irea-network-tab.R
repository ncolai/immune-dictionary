make_matrix <- function(df,rownames = NULL){
  # make matrix from data frame
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}

observeEvent(input$submit_cytokine_network, {
  runjs("$('.busy').show();")
  hideElement("plot_C")
  hideElement("download_C")
  hideElement("table_C")
  hideElement("download_table_C")
  #toggle('network_genediff_cutoff')
}, priority = 1, ignoreInit = TRUE)

observeEvent(input$dropdown_btn3,{
  # hide/show optional parameters on click
  toggle('network_genediff_cutoff')
})

genesC <- observeEvent(input$submit_cytokine_network, {
  `%notin%` = Negate(`%in%`)
  # hide tab A results
  hideElement("plot")
  hideElement("download")
  hideElement("table")
  hideElement("download_table")
  
  # hide tab B results
  hideElement("plot_B")
  hideElement("download_B")
  hideElement("table_B")
  hideElement("download_table_B")
  hideElement("radio_btns_B")
  
  if ((!is.null(input$network_matrix_singlefile) | input$sample_network_matrix_singlefile)){
    if (input$sample_network_matrix_singlefile){
      # use sample network tumor matrix data file
      my_data <- read_excel('exampleFiles/network_test.xlsx')
    }
    else{
      # use uploaded file
      table_data <- input$network_matrix_singlefile$datapath
      my_data <- read_excel(table_data, .name_repair = ~ ifelse(nzchar(.x), .x, paste('Sample', seq_along(.x) - 1)))  # rename columns without names
    }
    cat('calculating')
    data$input_profile_C <- make_matrix(select(my_data,-1), pull(my_data,1))      # make matrix from input
    
    cells_to_analyze = c() #vector of cells we want to analyze
    #if (input$sample_network_b_cell) {cells_to_analyze <- c(cells_to_analyze, "B_cell")}
    if (input$sample_network_t_cell_cd4) {cells_to_analyze <- c(cells_to_analyze, "T_cell_CD4")}
    if (input$sample_network_t_cell_cd8) {cells_to_analyze <- c(cells_to_analyze, "T_cell_CD8")}
    #if (input$sample_network_t_cell_gd) {cells_to_analyze <- c(cells_to_analyze, "T_cell_gd")}
    if (input$sample_network_treg) {cells_to_analyze <- c(cells_to_analyze, "Treg")}
    if (input$sample_network_nk_cell) {cells_to_analyze <- c(cells_to_analyze, "NK_cell")}
    if (input$sample_network_pdc) {cells_to_analyze <- c(cells_to_analyze, "pDC")}
    if (input$sample_network_cdc1) {cells_to_analyze <- c(cells_to_analyze, "cDC1")}
    if (input$sample_network_cdc2) {cells_to_analyze <- c(cells_to_analyze, "cDC2")}
    if (input$sample_network_migdc) {cells_to_analyze <- c(cells_to_analyze, "MigDC")}
    #if (input$sample_network_langerhans) {cells_to_analyze <- c(cells_to_analyze, "Langerhans")}
    if (input$sample_network_macrophage) {cells_to_analyze <- c(cells_to_analyze, "Macrophage")}
    #if (input$sample_network_monocyte) {cells_to_analyze <- c(cells_to_analyze, "Monocyte")}
    if (input$sample_network_neutrophil) {cells_to_analyze <- c(cells_to_analyze, "Neutrophil")}
    
    
    # check if file format is correct
    # make sure each column has a title and only numbers
    if (!isTRUE(all(is.numeric(data$input_profile_C) & -1000000 <= data$input_profile_C & data$input_profile_C <= 100000))){
      data$input_profile_C <- NULL
      showNotification(HTML(paste0('Your data file format is incorrect.', '<br>', 'Please make sure the first column is genes, 
                                   and the following columns only contain numbers (in the range [-1 to 1]) pertaining to the sample.')), duration = NULL, type = 'error')
    }
    else{
      # check if genes are valid
      # Parameter removed: ", species = tolower(input$speciesInput)"
      is_valid <- valid_genes(rownames(data$input_profile_C))
      
      if (identical(is_valid, 'invalid')){
        showNotification('Please submit valid genes in the first column.', duration = NULL, type = 'error')
        data$input_profile_C <- NULL
        data$table_tabC = NULL
        data$irea_plot_C = NULL
      }
      else {
        # if there are invalid genes, show error message
        if (!identical(is_valid, 'valid')){
          #note that there is a significant discrepancy because tumor data is from a different dataset
          #showNotification(HTML(paste0('The following genes were not found:', '<br>', paste(is_valid, collapse = ', '), '<br>', 'The calculation will proceed without these genes.', '<br><br>', 'Example mouse genes: Isg15, Irf7, Ncr1', '<br>', 'Example human genes: ISG15, IRF7, NCR1')), duration = 10, type = 'error')
        }
        cat('Calculating Cytokine Network Analysis...\n')
        tryCatch({
          # perform GetEnrichmentScoreProjection and IreaCompassPlot calculations
          cat("Input Cytokine Network\n")
          
          # To get cell-cell interaction
          library(openxlsx)
          df_irea_all = IreaAll(data$input_profile_C, species = "mouse", threshold_receptor = 0.05, threshold_ligand = 0.05, celltypes=cells_to_analyze)
          #save(df_irea_all, file = "df_irea_all.RData")
          #load("df_irea_all.RData")
          df_irea_network = IreaNetwork(df_irea_all, require_receptor_expression = TRUE)
          df_irea_network_pos = subset(df_irea_network, ES > 0)
          if (nrow(df_irea_network_pos) == 0) { #not enough connection data found, so throwing error
            stop("No significant connection data found for selected combination, please try again.")
          }
          
          data$table_tabC <- df_irea_network_pos
          
          cat("Got data table\n")
          # df_irea_pd1 = subset(data$table_tabB, Sample == colnames(data$input_profile_C)[1])
          data$irea_plot_C_type <- "Network"
          
          # OR
          # data$irea_plot_C <- IreaCompassPlot(data$table_tabC, color_by = "pval")
          
          cat("Plot IREA\n")
        },
        error = function(e){
          showNotification(paste0("Error: ", e), duration = NULL, type = 'error')
          #print(e)
          data$input_profile_C = NULL
          data$table_tabC = NULL
          data$irea_plot_C = NULL
        })
      }
    }
  }
  else{
    showNotification('Please upload gene matrix file', type = 'error')
  }
  runjs("$('.busy').hide();")
  showElement("plot_C")
  showElement("download_C")
  showElement("table_C")
  showElement("download_table_C")
}, ignoreInit = TRUE)

output$table_C <- renderDataTable({
  # if table data has been calculated, display table
  # req(data$table_tabC)
  # data$table_tabC
  
  # TODO: in the future, can select one of the samples
  req(data$table_tabC)
  # data$table_tabC_subset <- subset(data$table_tabC, Sample == input$rb)
  # data$table_tabC_subset
  data$table_tabC_subset <- data$table_tabC
})


plotInput_C <- function(){
  # if table data has been calculated, display plot
  # req(data$irea_plot_C)
  # data$irea_plot_C
  
  req(data$table_tabC)
  
  # if (!is.null(input$rb)){
  #   df_irea_pd1 = subset(data$table_tabC, Sample == input$rb)
    
    if (data$irea_plot_C_type == "Network") {
      data$irea_plot_C <- IreaNetworkCircosAll(data$table_tabC)
    }
    else if (data$irea_plot_C_type == "Network Individual") {
      data$irea_plot_C <- IreaNetworkCircosIndividual(data$table_tabC, n_cols = 4) 
    }
  # }
  
  data$irea_plot_C
}

output$plot_C <- renderPlot({
  plotInput_C()
}, height = 500)


output$download_network_matrix <- downloadHandler(
  # download example gene_matrix file
  filename = function(){"ex_network.xlsx"},
  content <- function(file) {
    file.copy("exampleFiles/network_test.xlsx", file)
  },
  contentType = "application/xlsx"
)

# observeEvent(input$network_matrix_singlefile, {
#   updateCheckboxInput(session, 'sample_network_matrix_singlefile', value = 'FALSE')
# })

observeEvent(input$sample_network_matrix_singlefile, {
  # session$sendCustomMessage("upload_txt", "SOME OTHER TEXT")
  if (input$sample_network_matrix_singlefile == TRUE){
    hideElement('network_matrix_singlefile')
  }
  else{
    showElement('network_matrix_singlefile')
  }
})

# download button tab B
output$download_C <- renderUI({
  req(data$irea_plot_C)
  tagList(downloadButton('downloadPDF_C', 'Download as PDF'),
          downloadButton('downloadPNG_C', 'Download as PNG'),
          downloadButton('downloadJPG_C', 'Download as JPG'),
          HTML("<br><br><br>"))
})

output$download_table_C <- renderUI({
  req(data$table_tabC)
  tagList(
    downloadButton('downloadTable_C', "Download table"),
    #downloadButton('downloadTable_C_all', "Download all samples"),
    HTML("<br><br><br>")
  )
})

output$downloadPDF_C <- downloadHandler(
  filename = "IREA_plot.pdf",
  content = function(file) {
    cairo_pdf(file, 6, 6)
    print(plotInput_C())
    dev.off()
  }) 

output$downloadPNG_C <- downloadHandler(
  filename = "IREA_plot.png",
  content = function(file) {
    png(file, width     = 6,
        height    = 6,
        units     = "in",
        res       = 1200,
        pointsize = 8)
    print(plotInput_C())
    dev.off()
  })    

output$downloadJPG_C <- downloadHandler(
  filename = "IREA_plot.jpg",
  content = function(file) {
    jpeg(file, width     = 6,
         height    = 6,
         units     = "in",
         res       = 1200,
         pointsize = 8)
    print(plotInput_C())
    dev.off()
  })    


# output$downloadTable_C <- downloadHandler(
#   filename = function(){"df_irea.csv"},
#   content = function(fname){
#     write.csv(data$table_tabB_subset, fname)
#   }
# )

output$downloadTable_C <- downloadHandler(
  filename = function(){"IREA_output.csv"},
  content = function(fname){
    write.csv(data$table_tabC, fname)
  }
)
