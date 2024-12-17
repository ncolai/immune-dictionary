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
  suppressWarnings({
    start_time <- Sys.time()
    
    if (input$network_type_cytokines == 'individual') { # & length(input$network_inputCytokines) > 0) {
      if (!is.null(input$network_matrix_singlefile) | input$sample_network_matrix_singlefile) {
        if (input$sample_network_matrix_singlefile) {
          # use sample network tumor matrix data file
          my_data <- read_excel('exampleFiles/network_test.xlsx')
        } else {
          # use uploaded file
          table_data <- input$network_matrix_singlefile$datapath
          my_data <- read_excel(table_data, .name_repair = ~ ifelse(nzchar(.x), .x, paste('Sample', seq_along(.x) - 1)))  # rename columns without names
        }
        
        cat('calculating')
        #data$my_data <- subset(my_data, "Cytokine" %in% input$network_inputCytokines)
        data$input_profile_C <- make_matrix(select(my_data, -1), pull(my_data, 1)) # make matrix from input
        
        cells_to_analyze <- c() #vector of cells we want to analyze
        if (input$sample_network_t_cell_cd4) {cells_to_analyze <- c(cells_to_analyze, "T_cell_CD4")}
        if (input$sample_network_t_cell_cd8) {cells_to_analyze <- c(cells_to_analyze, "T_cell_CD8")}
        if (input$sample_network_treg) {cells_to_analyze <- c(cells_to_analyze, "Treg")}
        if (input$sample_network_nk_cell) {cells_to_analyze <- c(cells_to_analyze, "NK_cell")}
        if (input$sample_network_pdc) {cells_to_analyze <- c(cells_to_analyze, "pDC")}
        if (input$sample_network_cdc1) {cells_to_analyze <- c(cells_to_analyze, "cDC1")}
        if (input$sample_network_cdc2) {cells_to_analyze <- c(cells_to_analyze, "cDC2")}
        if (input$sample_network_migdc) {cells_to_analyze <- c(cells_to_analyze, "MigDC")}
        if (input$sample_network_macrophage) {cells_to_analyze <- c(cells_to_analyze, "Macrophage")}
        if (input$sample_network_neutrophil) {cells_to_analyze <- c(cells_to_analyze, "Neutrophil")}
        
        # Check if file format is correct
        if (!isTRUE(all(is.numeric(data$input_profile_C) & -1000000 <= data$input_profile_C & data$input_profile_C <= 100000))) {
          data$input_profile_C <- NULL
          showNotification(HTML('Your data file format is incorrect.<br>Please make sure the first column is genes, and the following columns only contain numbers (in the range [-1 to 1]) pertaining to the sample.'), duration = NULL, type = 'error')
        } else {
          is_valid <- valid_genes(rownames(data$input_profile_C))
          if (identical(is_valid, 'invalid')) {
            showNotification('Please submit valid genes in the first column.', duration = NULL, type = 'error')
            data$input_profile_C <- NULL
            data$table_tabC = NULL
            data$irea_plot_C = NULL
          } else {
            if (!identical(is_valid, 'valid')) {
              # handle invalid genes message if needed
            }
            cat('Calculating Cytokine Network Analysis...\n')
            tryCatch({
              cat("Input Cytokine Network\n")
              #df_irea_all <- IreaAll(data$input_profile_C, species = "mouse", threshold_receptor = 0.05, threshold_ligand = 0.05, celltypes = cells_to_analyze)
              df_irea_all <- IreaAll(data$input_profile_C, species = input$speciesInputC, threshold_receptor = 0.05, threshold_ligand = 0.05, celltypes = cells_to_analyze)
              df_irea_network <- IreaNetwork(df_irea_all, require_receptor_expression = TRUE)
              df_irea_network_pos <- subset(df_irea_network, ES > 0)
              if (nrow(df_irea_network_pos) == 0) {
                stop("No significant connection data found for selected combination, please try again.")
              }
              data$table_tabC <- df_irea_network_pos
              cat("Got data table\n")
            },
            error = function(e) {
              showNotification(paste0("Error: ", e), duration = NULL, type = 'error')
              data$input_profile_C = NULL
              data$table_tabC = NULL
              data$irea_plot_C = NULL
            })
          }
        }
      } else {
        showNotification('Please upload gene matrix file', type = 'error')
      }
      runjs("$('.busy').hide();")
      showElement("plot_C")
      showElement("download_C")
      showElement("table_C")
      showElement("download_table_C")
    } else if (input$network_type_cytokines == 'inp' & length(input$network_inputCytokines) <= 0) {
      showNotification('Please enter cytokines.', type = 'error')
    } else if (input$network_type_cytokines == 'all') {
      if (!is.null(input$network_matrix_singlefile) | input$sample_network_matrix_singlefile) {
        if (input$sample_network_matrix_singlefile) {
          # use sample network tumor matrix data file
          my_data <- read_excel('exampleFiles/network_test.xlsx')
        } else {
          # use uploaded file
          table_data <- input$network_matrix_singlefile$datapath
          my_data <- read_excel(table_data, .name_repair = ~ ifelse(nzchar(.x), .x, paste('Sample', seq_along(.x) - 1)))  # rename columns without names
        }
        
        cat('calculating')
        data$input_profile_C <- make_matrix(select(my_data, -1), pull(my_data, 1)) # make matrix from input
        
        cells_to_analyze <- c() #vector of cells we want to analyze
        if (input$sample_network_t_cell_cd4) {cells_to_analyze <- c(cells_to_analyze, "T_cell_CD4")}
        if (input$sample_network_t_cell_cd8) {cells_to_analyze <- c(cells_to_analyze, "T_cell_CD8")}
        if (input$sample_network_treg) {cells_to_analyze <- c(cells_to_analyze, "Treg")}
        if (input$sample_network_nk_cell) {cells_to_analyze <- c(cells_to_analyze, "NK_cell")}
        if (input$sample_network_pdc) {cells_to_analyze <- c(cells_to_analyze, "pDC")}
        if (input$sample_network_cdc1) {cells_to_analyze <- c(cells_to_analyze, "cDC1")}
        if (input$sample_network_cdc2) {cells_to_analyze <- c(cells_to_analyze, "cDC2")}
        if (input$sample_network_migdc) {cells_to_analyze <- c(cells_to_analyze, "MigDC")}
        if (input$sample_network_macrophage) {cells_to_analyze <- c(cells_to_analyze, "Macrophage")}
        if (input$sample_network_neutrophil) {cells_to_analyze <- c(cells_to_analyze, "Neutrophil")}
        
        # Check if file format is correct
        if (!isTRUE(all(is.numeric(data$input_profile_C) & -1000000 <= data$input_profile_C & data$input_profile_C <= 100000))) {
          data$input_profile_C <- NULL
          showNotification(HTML('Your data file format is incorrect.<br>Please make sure the first column is genes, and the following columns only contain numbers (in the range [-1 to 1]) pertaining to the sample.'), duration = NULL, type = 'error')
        } else {
          is_valid <- valid_genes(rownames(data$input_profile_C))
          if (identical(is_valid, 'invalid')) {
            showNotification('Please submit valid genes in the first column.', duration = NULL, type = 'error')
            data$input_profile_C <- NULL
            data$table_tabC = NULL
            data$irea_plot_C = NULL
          } else {
            if (!identical(is_valid, 'valid')) {
              # handle invalid genes message if needed
            }
            cat('Calculating Cytokine Network Analysis...\n')
            tryCatch({
              cat("Input Cytokine Network\n")
              df_irea_all <- IreaAll(data$input_profile_C, species = input$speciesInputC, threshold_receptor = 0.05, threshold_ligand = 0.05, celltypes = cells_to_analyze)
              df_irea_network <- IreaNetwork(df_irea_all, require_receptor_expression = TRUE)
              df_irea_network_pos <- subset(df_irea_network, ES > 0)
              if (nrow(df_irea_network_pos) == 0) {
                stop("No significant connection data found for selected combination, please try again.")
              }
              data$table_tabC <- df_irea_network_pos
              cat("Got data table\n")
            },
            error = function(e) {
              showNotification(paste0("Error: ", e), duration = NULL, type = 'error')
              data$input_profile_C = NULL
              data$table_tabC = NULL
              data$irea_plot_C = NULL
            })
          }
        }
      } else {
        showNotification('Please upload gene matrix file', type = 'error')
      }
      runjs("$('.busy').hide();")
      showElement("plot_C")
      showElement("download_C")
      showElement("table_C")
      showElement("download_table_C")
    }
    
    end_time <- Sys.time()
    elapsed_time <- end_time - start_time
    print(elapsed_time)
  })
})

output$table_C <- renderDataTable({
  req(data$table_tabC)
  data$table_tabC_subset <- data$table_tabC
})


plotInput_C <- function(){
  req(data$table_tabC)
    if (input$network_type_cytokines == "all") {
      data$irea_plot_C <- IreaNetworkCircosAll(data$table_tabC)
    }
    else if (input$network_type_cytokines == "individual") {
      data$irea_plot_C <- IreaNetworkCircosIndividual(data$table_tabC, n_cols = 4)#, input$network_inputCytokines) 
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

output$downloadTable_C <- downloadHandler(
  filename = function(){"IREA_output.csv"},
  content = function(fname){
    write.csv(data$table_tabC, fname)
  }
)

# output$download_network_all <- renderUI({
#   req(data$network_plot_all)
#   tagList(downloadButton('download_networks_all_csv', 'Download Network Plot as CSV'),
#           downloadButton('download_networks_all_pdf', 'Download Network Plot as PDF'),
#           HTML("<br><br><br>"))
# })
# 
# output$download_network_inp <- renderUI({
#   req(data$network_plot_input)
#   tagList(downloadButton('download_networks_inp_csv', 'Download Network Plot as CSV'),
#           downloadButton('download_networks_inp_pdf', 'Download Network Plot as PDF'),
#           HTML("<br><br><br>"))
# })
# 
# output$download_networks_all_csv <- downloadHandler(
#   filename = function(){"network_plot_all.csv"},
#   content = function(fname){
#     write.csv(data$network_all_data, fname)
#   }
# )
# 
# output$download_networks_inp_csv <- downloadHandler(
#   filename = function(){"network_plot_inp.csv"},
#   content = function(fname){
#     write.csv(data$network_input_data, fname)
#   }
# )
# 
# output$download_networks_all_pdf <- downloadHandler(
#   filename = function(){paste("network_plot_all ",Sys.Date(),".pdf",sep="")},
#   content = function(file){
#     saveWidget(data$violin_plot_all, 'temp.html', selfcontained = FALSE)
#     webshot2::webshot('temp.html', file = file, zoom = 5)
#     unlink('temp.html')
#   }
# )
# 
# output$download_networks_inp_pdf <- downloadHandler(
#   filename = function(){paste("network_plot_inp ",Sys.Date(),".pdf",sep="")},
#   content = function(file){
#     saveWidget(data$violin_plot_inp, 'temp.html', selfcontained = FALSE)
#     webshot2::webshot('temp.html', file = file, zoom = 5)
#     unlink('temp.html')
#   }
# )
