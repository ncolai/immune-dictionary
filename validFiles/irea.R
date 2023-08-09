library(plotly)
library(tidyverse)
library(shinyBS)
library(shinyWidgets)
library(readxl)
library(shinyjs)
library(tidyverse)
source('sourceFiles/irea-function.R')
source('sourceFiles/irea-visualization.R')
source('sourceFiles/global.R')
source('sourceFiles/irea-gene-list-tab.R', local = TRUE)
source('sourceFiles/irea-gene-matrix-tab.R', local = TRUE)

cellList <- c(" ", readLines("sourceFiles/lig_seurat_data.txt"))

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
  titlePanel("IREA"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(id = "tabs",
        tabPanel("Gene List", value = "A", 
                 h3(),
                 bsTooltip("info_btn", "A list of genes you would like to investigate, such as differentially expressed genes obtained from DESeq2 or Seurat FindMarkers", 
                           placement = "right", trigger = "hover", options = NULL),
                 textAreaInput(inputId="inputGene", label = tags$div(class = "header", tags$div("Input Genes", tags$a(style = "margin-left:0.2em; text-decoration: none;", id = "info_btn", 'ⓘ'))), rows = 15, placeholder = 'Gene1\nGene2\n...', resize = "vertical"),
                 actionButton("sample", "Sample Genes", class = "btn-link"),
                 actionButton("clear", "Clear Genes", class = "btn-link"),
                 div(class = 'gene_file_text', style=
                       'margin: auto;
                     text-align: center;',
                                 '- or -'),
                 fileInput("gene_file", "Upload Gene File",
                           multiple = FALSE,
                           accept = c(".txt", ".xlsx", ".xls")),
                 #HTML("<br>"),
                 selectInput("inputCell", "Choose Cell Type", cellList),
                 div(style='overflow-x: auto;',
                     radioGroupButtons(inputId = "inputMethod", label = "Method", choices = c("Score", "Hypergeometric"), justified = TRUE)),
                 HTML("<br>"),
                 actionButton("submit", "Submit", class = "btn-block"),
                ),
        tabPanel("Gene Matrix", value = "B",
                 h3(),
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
                 selectInput("inputCell_tabB", "Choose Cell Type", cellList),
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
                 hidden(numericInput("genediff_cutoff", "Gene Diff Cutoff", 0.25, min = 0, max = 1, step = 0.05)),
                 HTML("<br>"),
                 actionButton("submit_tabB", "Submit", class = "btn-block"),
                )
      ),
      width = 3
    ),
    mainPanel(
      tags$div(class = "busy", 
               tags$img(src = "https://media.giphy.com/media/sSgvbe1m3n93G/giphy.gif")),
      # tab A
      plotOutput('plot', height = '500px'),
      uiOutput("download"),
      dataTableOutput("table"),
      uiOutput("download_table"),
      
      # tab B
      uiOutput('radio_btns_B'),
      plotOutput('plot_B', height = '500px'),
      uiOutput("download_B"),
      dataTableOutput('table_B'),
      uiOutput("download_table_B")
      #downloadButton('downloadTable_B_all', "Download all samples")
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

data <- reactiveValues(table = NULL, compass_plot = NULL, table_tabB = NULL, compass_plot_B = NULL, input_profile = NULL, table_tabB_subset = NULL)

"%ni%" <- Negate("%in%")
valid_genes <- function(gene_list) {
  # check if list of genes is valid (in lig_seurat)
  not_working_genes <- vector()
  final_genes_list <- vector()
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
}

split_input <- function(input_list){
  # split input into vector output
  split_vec <- strsplit(input_list, "\\s*,\\s*|\\s+")
  return (split_vec[[1]])
}