library(plotly)
library(tidyverse)
library(shinyBS)
library(shinyWidgets)
library(readxl)
library(shinyjs)
library(tidyverse)

cytokineList <- c(readLines("sourceFiles/cytokine_list.txt"), "")
newCytokineList <- subset(read.xlsx("dataFiles/irea_cytokine_list.xlsx"), select = "Cytokine_DisplayName")


source('sourceFiles/irea-function.R')
source('sourceFiles/irea-visualization.R')
source('sourceFiles/irea-receptor.R')
source('sourceFiles/irea-network.R')
source('sourceFiles/global.R')
source('sourceFiles/irea-gene-list-tab.R', local = TRUE)
source('sourceFiles/irea-gene-matrix-tab.R', local = TRUE)
source('sourceFiles/irea-network-tab.R', local = TRUE)


cellListIreaList <- c(" ", readLines("sourceFiles/lig_seurat_listcells.txt"))
cellListIreaMatrix <- c(" ", readLines("sourceFiles/lig_seurat_matrixcells.txt"))

css <- "
.busy { 
  position: fixed;
  z-index: 1000;
  top: 50%;
  left: 50%;
  margin-top: -100px;
  margin-left: -50px;
  display: none;
  text-align: center;
  padding-top: 20px;
  padding-left: 30px;
  padding-bottom: 40px;
  padding-right: 30px;
  border-radius: 5px;
}

.shiny-notification {
          position:fixed;
          top: calc(22%);
          left: calc(25%);
        }
"

# jscode_upload_txt <- " Shiny.addCustomMessageHandler('upload_txt', function(txt) {   
# $('input:text.form-control').attr('placeholder','Some New Text'); }); "

output$pageStub <- renderUI(fluidPage(
  #' tags$head(
  #'   # Note the wrapping of the string in HTML()
  #'   tags$style(HTML("
  #'     @import url('https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css');
  #'     "))
  #' ),
  #tags$style(HTML('<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css">')),
  useShinyjs(),
  tags$head(tags$style(css)),
  titlePanel("Immune Response Enrichment Analysis"),
  h4("Analyze your own gene expression data"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(id = "tabs",
                  tabPanel("Gene List", value = "A", 
                           h3(),
                           bsTooltip("info_btn", HTML(paste0("A list of genes you would like to investigate, such as differentially expressed genes obtained from DESeq2 or Seurat FindMarkers.")), 
                                     placement = "right", trigger = "hover", options = NULL),
                           textAreaInput(inputId="inputGene", label = tags$div(class = "header", tags$div("Input Genes", tags$a(style = "margin-left:0.2em; text-decoration: none;", id = "info_btn", 'ⓘ'))), rows = 10, placeholder = 'Gene1\nGene2\n...', resize = "vertical"),
                           actionButton("sample", "Sample Genes", class = "btn-link"),
                           actionButton("clear", "Clear Genes", class = "btn-link"),
                           #div(class = 'gene_file_text', style=
                           #       'margin: auto;
                           #text-align: center;',
                           #   '- or -'),
                           #fileInput("gene_file", "Upload Gene File",
                           #           multiple = FALSE,
                           #           accept = c(".txt", ".xlsx", ".xls")),
                           HTML("<br>"),
                           HTML("<br>"),
                           selectInput("inputCell", "Choose Cell Type", cellListIreaList),
                           div(style='overflow-x: auto;',
                               radioGroupButtons(inputId = "speciesInput", label = "Species", choices = c("Mouse", "Human"), justified = TRUE)),
                           HTML('<hr>'),
                           actionButton("dropdown_btn", label = ' Optional Parameters', icon = icon("chevron-down"),
                                        class = "btn-link", style = 'text-decoration: none;
                  padding: 0;
                  margin: 0;
                  color: black;
                  font-weight: bold;
                  font-size: 1em;
                  white-space: normal;
                  width: 100%;
                  text-align: left;'),
                           HTML('<br><br>'),
                           hidden(radioGroupButtons(inputId = "inputMethod", label = "Method", choices = c("Score", "Hypergeometric"), justified = TRUE)),
                           #  hidden(radioGroupButtons(inputId = "speciesInputB", label = "Species", choices = c("Mouse", "Human"), justified = TRUE)),
                           HTML("<hr>"),
                           h3("Submit"),
                           actionButton("submit_radar_list", "Analyze Cell Polarization", class = "btn-block"),
                           HTML("<br>"),
                           actionButton("submit_compass_list", "Analyze Cytokine Response", class = "btn-block"),
                           
                           HTML('<hr>'), 
                           p("To interpret cell polarization analysis output, please refer to the \"Reference Cards\" tab above for detailed info of each polarization state.")
                           
                           # HTML("<br>"),
                           # actionButton("submit_network_list", "To be announced!", class = "btn-block")
                  ),
                  tabPanel("Gene Matrix", value = "B",
                           h3(),
                           HTML('<br>'),     
                           # tags$script(HTML(jscode_upload_txt)),
                           div(style=
                                 'padding: 0;
                        margin: 0;
                        color: black;
                        font-weight: bold;
                        font-size: 1em;
                        white-space: normal;
                        width: 100%;
                        text-align: left;
                        margin-bottom: 5px;',
                               HTML('Upload Gene Matrix File')
                           ),
                           fileInput("matrix_file", label = NULL,
                                     multiple = FALSE,
                                     accept = c(".txt", ".xlsx", ".xls")),
                           materialSwitch(
                             inputId = "sample_matrix",
                             label = "Use Example File", 
                             status = "primary",
                             right = TRUE
                           ),
                           # checkboxInput(inputId = "sample_matrix", label = "Use Example File"),
                           
                           # actionButton("sample_matrix", "Use Example File", class = "btn-link"),
                           div(style='overflow-x: auto;',
                               downloadButton("download_matrix", "Download Example File", class = "btn-link")),
                           HTML('<br>'),
                           h3(),
                           bsTooltip("info_btn2", "Please use NK Cell for example files", 
                                     placement = "right", trigger = "hover", options = NULL),
                           selectInput("inputCell_tabB", label = tags$div(class = "header", tags$div("Choose Cell Type", tags$a(style = "margin-left:0.2em; text-decoration: none;", id = "info_btn2", 'ⓘ'))), cellListIreaMatrix),
                           HTML('<hr>'),
                           actionButton("dropdown_btn2", label = ' Optional Parameters', icon = icon("chevron-down"),
                                        class = "btn-link", style = 'text-decoration: none;
                  padding: 0;
                  margin: 0;
                  color: black;
                  font-weight: bold;
                  font-size: 1em;
                  white-space: normal;
                  width: 100%;
                  text-align: left;'),
                           HTML('<br><br>'),
                           hidden(numericInput("genediff_cutoff", "Gene Diff Cutoff", 0.25, min = 0, max = 1, step = 0.05)),
                           #  hidden(radioGroupButtons(inputId = "speciesInputB", label = "Species", choices = c("Mouse", "Human"), justified = TRUE)),
                           HTML("<hr>"),
                           HTML("<br>"),
                           h3("Submit"),
                           actionButton("submit_radar_matrix", "Analyze Cell Polarization", class = "btn-block"),
                           HTML("<br>"),
                           actionButton("submit_compass_matrix", "Analyze Cytokine Response", class = "btn-block"),
                           HTML('<hr>'), 
                           p("To interpret cell polarization analysis output, please refer to the \"Reference Cards\" tab above for detailed info of each polarization state."),
                           HTML('<hr>'), 
                           p("It takes 5-10min to run the gene matrix input analysis. Simply leave the browser on and come back later to see the results.")
                           
                           # HTML("<br>"),
                           # actionButton("submit_network_matrix", "To be announced!", class = "btn-block")
                  )
                  , 
                  tabPanel("Network", value = "C",
                           h3(),
                           HTML('<br>'),     
                           
                           # tags$script(HTML(jscode_upload_txt)),
                           div(style=
                                 'padding: 0;
                        margin: 0;
                        color: black;
                        font-weight: bold;
                        font-size: 1em;
                        white-space: normal;
                        width: 100%;
                        text-align: left;
                        margin-bottom: 5px;',
                               HTML('Upload an Excel file containing the gene information for Cytokine Network Analysis. See example for how to format excel file.')
                           ),
                           
                           fileInput("network_matrix_singlefile", label = NULL,
                                     multiple = FALSE,
                                     accept = c(".txt", ".xlsx", ".xls")),
                           materialSwitch(
                             inputId = "sample_network_matrix_singlefile",
                             label = "Use Example File", 
                             status = "primary",
                             right = TRUE
                           ),

                           div(style='overflow-x: auto;',
                               downloadButton("download_network_matrix", "Download Example File", class = "btn-link")),
                           HTML('<hr>'),
                           div(style=
                                 'padding: 0;
                        margin: 0;
                        color: black;
                        font-weight: bold;
                        font-size: 1em;
                        white-space: normal;
                        width: 100%;
                        text-align: left;
                        margin-bottom: 5px;',
                               HTML('Choose 2 or more cells to analyze'),
                           ),
                           fluidRow(
                             # First column
                             column(6,
                                    #h4("B cell"),
                                    #materialSwitch(
                                    #  inputId = "sample_network_b_cell",
                                    #  label = "Analyze B Cell", 
                                    #  status = "primary",
                                    #  right = TRUE
                                    #),
                                    
                                    h4("CD4+ T cell"),
                                    materialSwitch(
                                      inputId = "sample_network_t_cell_cd4",
                                      label = "Analyze CD4+ T Cell", 
                                      status = "primary",
                                      right = TRUE
                                    ),
                                    
                                    h4("CD8+ T cell"),
                                    materialSwitch(
                                      inputId = "sample_network_t_cell_cd8",
                                      label = "Analyze CD8+ T Cell", 
                                      status = "primary",
                                      right = TRUE
                                    ),
                                    
                                    #h4("γδ T cell"),
                                    #materialSwitch(
                                    #  inputId = "sample_network_t_cell_gd",
                                    #  label = "Analyze gamma-delta T Cell", 
                                    #  status = "primary",
                                    #  right = TRUE
                                    #),
                                    
                                    h4("Treg"),
                                    materialSwitch(
                                      inputId = "sample_network_treg",
                                      label = "Analyze Treg Cell", 
                                      status = "primary",
                                      right = TRUE
                                    ),
                                    
                                    h4("NK cell"),
                                    materialSwitch(
                                      inputId = "sample_network_nk_cell",
                                      label = "Analyze NK Cell", 
                                      status = "primary",
                                      right = TRUE
                                    ),
                                    
                                    h4("pDC"),
                                    materialSwitch(
                                      inputId = "sample_network_pdc",
                                      label = "Analyze pDC Cell", 
                                      status = "primary",
                                      right = TRUE
                                    ),
                                    HTML('<hr>')
                             ), 
                             # Second column
                             column(6,
                                    h4("cDC1"),
                                    materialSwitch(
                                      inputId = "sample_network_cdc1",
                                      label = "Analyze cDC1 Cell", 
                                      status = "primary",
                                      right = TRUE
                                    ),
                                    h4("cDC2"),
                                    materialSwitch(
                                      inputId = "sample_network_cdc2",
                                      label = "Analyze cDC2 Cell", 
                                      status = "primary",
                                      right = TRUE
                                    ),
                                    h4("MigDC"),
                                    materialSwitch(
                                      inputId = "sample_network_migdc",
                                      label = "Analyze MigDC Cell", 
                                      status = "primary",
                                      right = TRUE
                                    ),
                                    #h4("Langerhans"),
                                    #materialSwitch(
                                    #  inputId = "sample_network_langerhans",
                                    #  label = "Analyze Langerhans Cell", 
                                    #  status = "primary",
                                    #  right = TRUE
                                    #),
                                    h4("Macrophage"),
                                    materialSwitch(
                                      inputId = "sample_network_macrophage",
                                      label = "Analyze Macrophage Cell", 
                                      status = "primary",
                                      right = TRUE
                                    ),
                                    #h4("Monocyte"),
                                    #materialSwitch(
                                    #  inputId = "sample_network_monocyte",
                                    #  label = "Analyze Monocyte Cell", 
                                    #  status = "primary",
                                    #  right = TRUE
                                    #),
                                    h4("Neutrophil"),
                                    materialSwitch(
                                      inputId = "sample_network_neutrophil",
                                      label = "Analyze Neutrophil Cell", 
                                      status = "primary",
                                      right = TRUE
                                    ),
                             )
                           ),

                           # fluidRow(
                           #   # First column
                           #   column(6,
                           #          h4("B cell"),
                           #          fileInput("network_matrix_file_B_cell", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #          
                           #          h4("CD4+ T cell"),
                           #          fileInput("network_matrix_file_T_cell_CD4", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #          
                           #          h4("CD8+ T cell"),
                           #          fileInput("network_matrix_file_T_cell_CD8", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #          
                           #          h4("γδ T cell"),
                           #          fileInput("network_matrix_file_T_cell_gd", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #          
                           #          h4("Treg"),
                           #          fileInput("network_matrix_file_Treg", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #          
                           #          h4("NK cell"),
                           #          fileInput("network_matrix_file_NK_cell", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #          
                           #          h4("pDC"),
                           #          fileInput("network_matrix_file_pDC", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #          HTML('<hr>')
                           #   ),
                           #   
                           #   # Second column
                           #   column(6,
                           #          h4("cDC1"),
                           #          fileInput("network_matrix_file_cDC1", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #          h4("cDC2"),
                           #          fileInput("network_matrix_file_cDC2", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #          h4("MigDC"),
                           #          fileInput("network_matrix_file_MigDC", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #          h4("Langerhans"),
                           #          fileInput("network_matrix_file_Langerhans", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #          h4("Macrophage"),
                           #          fileInput("network_matrix_file_Macrophage", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #          h4("Monocyte"),
                           #          fileInput("network_matrix_file_Monocyte", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #          h4("Neutrophil"),
                           #          fileInput("network_matrix_file_Neutrophil", label = NULL,
                           #                    multiple = FALSE,
                           #                    accept = c(".txt", ".xlsx", ".xls"),
                           #                    placeholder = "Select file"),
                           #   )
                           # ),
                           actionButton("dropdown_btn3", label = ' Optional Parameters', icon = icon("chevron-down"),
                                        class = "btn-link", style = 'text-decoration: none;
                  padding: 0;
                  margin: 0;
                  color: black;
                  font-weight: bold;
                  font-size: 1em;
                  white-space: normal;
                  width: 100%;
                  text-align: left;'),
                           HTML('<br><br>'),
                           hidden(numericInput("network_genediff_cutoff", "Gene Diff Cutoff", 0.25, min = 0, max = 1, step = 0.05)),
                           #radioGroupButtons(inputId = "network_type_cytokines", 
                                             #label = "Choose Cytokine Input", 
                                             #choices = c("All", "Input"), justified = TRUE),
                           radioButtons(inputId="network_type_cytokines", label="Choose Cytokine Input",
                                        choices=c("All cytokines" = "all", "Individual cytokines" = "individual")),
                           #selectInput("network_inputCytokines", "Choose Cytokine Input", 
                                       #choices = cytokineList,
                                       #multiple = TRUE),
                           h3("Submit"),
                           actionButton("submit_cytokine_network", "Analyze Cytokine Network", class = "btn-block"),
                           HTML('<hr>'), 
                           p("It takes ~10min per cell type. Simply leave the browser on and come back later to see the results.")
                  )
                  
                  
      ),
      width = 4
    ),
    mainPanel(
      tags$div(class = "busy", 
               tags$img(src = "https://media.giphy.com/media/sSgvbe1m3n93G/giphy.gif")),
      # tab A
      plotOutput('plot', height = '500px'),
      uiOutput("download"),
      # dataTableOutput("table"),
      DT::DTOutput("table"),
      uiOutput("download_table"),
      
      # tab B
      uiOutput('radio_btns_B'),
      plotOutput('plot_B', height = '500px'),
      uiOutput("download_B"),
      # dataTableOutput('table_B'),
      DT::DTOutput('table_B'),
      uiOutput("download_table_B"),
      #downloadButton('downloadTable_B_all', "Download all samples")
      
      # tab C
      plotOutput('plot_C', height = '500px'),
      uiOutput("download_C"),
      # dataTableOutput('table_C'),
      DT::DTOutput('table_C'),
      uiOutput("download_table_C")
      #downloadButton('downloadTable_C_all', "Download all samples")
    ))))


observeEvent(input$sample,{
  updateTextAreaInput(session,
                      "inputGene",
                      value = "Ifitm3\nIsg15\nIfit3\nBst2\nSlfn5\nIsg20\nPhf11b\nZbp1\nRtp4"
  )
  return("Ifitm3\nIsg15\nIfit3\nBst2\nSlfn5\nIsg20\nPhf11b\nZbp1\nRtp4")
})

observeEvent(input$clear,{
  updateTextAreaInput(session, "inputGene", value = "")
  return("")
})

data <- reactiveValues(
  table = NULL, 
  irea_plot = NULL, 
  table_tabB = NULL, 
  irea_plot_B = NULL, 
  input_profile = NULL, 
  table_tabB_subset = NULL,
  irea_plot_type = NULL,
  input_profile_C = NULL, 
  table_tabC = NULL, 
  irea_plot_C = NULL,
  irea_plot_C_type = NULL)

"%ni%" <- Negate("%in%")
valid_genes <- function(gene_list, species = "mouse") {
  # check if list of genes is valid (in lig_seurat)
  not_working_genes <- vector()
  final_genes_list <- vector()
  
  if (species == "mouse") {
    for (gene in gene_list){
      if (gene %ni% valid_gene_list){
        not_working_genes <- c(not_working_genes, gene)
      }
      else {
        final_genes_list <- c(final_genes_list, gene)
      }
    }
    if (length(not_working_genes) == length(gene_list)){
      return('invalid')
    }
    else if (length(not_working_genes) > 0){
      return(not_working_genes)
    }
    else{
      return('valid')
    }
  } else {
    for (gene in gene_list){
      if (gene %ni% valid_gene_list_human){
        not_working_genes <- c(not_working_genes, gene)
      }
      else {
        final_genes_list <- c(final_genes_list, gene)
      }
    }
    if (length(not_working_genes) == length(gene_list)){
      return('invalid')
    }
    else if (length(not_working_genes) > 0){
      return(not_working_genes)
    }
    else{
      return('valid')
    }
  }
}

split_input <- function(input_list){
  # split input into vector output
  split_vec <- strsplit(input_list, "\\s*,\\s*|\\s+")
  return (split_vec[[1]])
}