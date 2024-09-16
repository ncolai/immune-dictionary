library(plotly)
library(tidyverse)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders)
library(reshape2)
library(ggplot2)
library(scales)
library(patchwork)
library(cowplot)
library(openxlsx)
library(fst)
#library(furrr)
#library(future)
library(webshot2)
library(htmlwidgets)

cellList <- c(readLines("sourceFiles/lig_seurat_data.txt"))
cytokineList <- c(readLines("sourceFiles/cytokine_list.txt"), "")
geneList <- c(readLines("sourceFiles/gene_list.txt"))
newCytokineList <- subset(read.xlsx("dataFiles/irea_cytokine_list.xlsx"), select = "Cytokine_DisplayName")
lig_seurat <- vector(mode="list", length=15)
names(lig_seurat) <- readLines("sourceFiles/lig_seurat_data.txt")
source('sourceFiles/data_browser_sidebar.R', local = TRUE)
source('sourceFiles/data_browser_autogenerate.R', local = TRUE)
source('sourceFiles/data_browser_heat_map.R', local = TRUE)
source('sourceFiles/data_browser_violin_plot.R', local = TRUE)
source('sourceFiles/data_browser_box_plot.R', local = TRUE)
source('sourceFiles/data_browser_umap_good.R', local = TRUE)

#if (!exists('example_plots')){
#  cat('Loading data\n')
#  #example_plots <- readRDS("dataFiles/example_plots.rds")
#  example_plots <- qread("dataFiles/example_plots.qs")
#  cat('Loaded data')
#}

css <- "
.shiny-notification {
          position:fixed;
          top: calc(22%);
          left: calc(25%);
          width: 250px;
        }
"
# observeEvent(input$submit, {
#   runjs("$('.busy').show();")
#   hide("plot")
# }, priority = 1, ignoreInit = TRUE)

data <- reactiveValues(plot_all1 = NULL, 
                       plot_input1 = NULL, 
                       plot_all2 = NULL, 
                       plot_input2 = NULL, 
                       violin_plot_all = NULL, 
                       violin_all_data = NULL,
                       violin_plot_input = NULL,
                       violin_input_data = NULL,
                       box_plot_all = NULL,
                       box_all_data = NULL,
                       box_plot_input = NULL,
                       box_input_data = NULL,
                       plot_res = NULL)

heatmap_data <- reactiveValues(selectedGenes = c("Il4i1", "St7", "Isg15"))
box_data <- reactiveValues(selectedGenes = c("Il4i1", "St7", "Isg15"))
violin_data <- reactiveValues(selectedGenes = c("Il4i1", "St7", "Isg15"))
umap_data <- reactiveValues(selectedGenes = NULL)

"%ni%" <- Negate("%in%")
valid_comp_cytokines <- function(gene_list, valid_genes) {
  not_working_genes <- vector()
  final_genes_list <- vector()
  for (gene in gene_list){
    if (gene %ni% valid_genes){
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
  split_vec <- strsplit(input_list, "\\s*,\\s*|\\s+")
  return (split_vec[[1]])
}

output$plot <- renderPlotly({
  if (is.null(input$tabs)) {
    plotly_empty()
  } else if (input$dataType == "Differential<br>Gene Expression") {
    showTab(inputId = "tabs", target = "A")
    showTab(inputId = "tabs", target = "B")
    hideTab(inputId = "tabs", target = "C")
    hideTab(inputId = "tabs", target = "D")
    hideTab(inputId = "tabs", target = "E")
    if (input$tabs == 'A'){
      if (input$featureInput1 == " ") return(NULL)
      else if (input$featureInput2 == " ") return(NULL)
      else{
        allsigs <- read.xlsx("dataFiles/SuppTable3_CytokineSignatures.xlsx")
        subset_data <- subset(allsigs, Celltype_Str == input$featureInput1)
        final_data <- subset(subset_data, Cytokine == input$featureInput2)
        gene_column <- final_data %>% pull(Gene)
        log_column <- final_data %>% pull(Avg_log2FC)
        
        if (length(gene_column) == 0) {
          plot_ly() %>%
            add_trace(type = "scatter", mode = "markers", x = c(0), y = c(0), marker = list(color = "rgba(0,0,0,0)"), showlegend = FALSE) %>%
            layout(
              annotations = list(
                text = "No DEG Data Available",
                showarrow = FALSE,
                x = 0.5,
                y = 0.5,
                xref = "paper",
                yref = "paper"
              ),
              xaxis = list(showline = FALSE, zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, range = c(NA, NA)),
              yaxis = list(showline = FALSE, zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE, range = c(NA, NA))
            )
        } else {
          y <- factor(gene_column)
          x <- log_column
          plot_height <- 500 + 15*nrow(final_data)
          
          fig <- plot_ly(
            y = reorder(y,x),
            x = x,
            name = "DEG Cytokine Expression",
            type = "bar",
            text = paste0("p-value: ",final_data %>% pull(FDR),sep=""),
            height = plot_height) %>%
            
            # x <- factor(gene_column)
            # y <- log_column
            # 
            # fig <- plot_ly(
            #   x = reorder(x,desc(y)),
            #   y = y,
            #   name = "DEG Cytokine Expression",
            #   type = "bar") %>%
            
            layout(title = 'DEG Cytokine Expression',
                   yaxis = list(title = 'Gene'),
                   xaxis = list(side = 'top',title = 'Magnitude of Differential Expression'),
                   margin = list(t=150, l=100))
          
          fig
        }
      }
    }
    else if (input$tabs == 'B'){
      updateSelectizeInput(session, 'inputGene', choices = geneList, server = TRUE, selected = heatmap_data$selectedGenes)
      if (input$inputMethod == "Cell Type"){
        shinyjs::show('cellType1')
        # show('inputCytokines')
        shinyjs::show('submit_celltype')
        shinyjs::show('type_cytokines')
        shinyjs::hide('cytokine1')
        shinyjs::hide('inputCellTypes')
        shinyjs::hide('submit_cytokines')
        shinyjs::hide('type_cellType')
        
        if (input$type_cytokines == 'all'){
          shinyjs::hide('inputCytokines')
          req(data$plot_all1)
          data$plot_all1
        }
        else{
          shinyjs::show('inputCytokines')
          req(data$plot_input1)
          data$plot_input1
        }
      }
      else {
        shinyjs::hide('cellType1')
        shinyjs::hide('inputCytokines')
        shinyjs::hide('submit_celltype')
        shinyjs::hide('type_cytokines')
        shinyjs::show('cytokine1')
        #show('inputCellTypes')
        shinyjs::show('submit_cytokines')
        shinyjs::show('type_cellType')
        
        if (input$type_cellType == 'all'){
          shinyjs::hide('inputCellTypes')
          req(data$plot_all2)
          data$plot_all2
        }
        else{
          shinyjs::show('inputCellTypes')
          req(data$plot_input2)
          data$plot_input2
        }
      }
    }
  }
  else if (input$dataType == "Gene<br>Expression") {
    hideTab(inputId = "tabs", target = "A")
    hideTab(inputId = "tabs", target = "B")
    showTab(inputId = "tabs", target = "C")
    showTab(inputId = "tabs", target = "D")
    showTab(inputId = "tabs", target = "E")
    if (input$tabs == 'C'){
      updateSelectizeInput(session, 'violinInputGene', choices = geneList, server = TRUE, selected = violin_data$selectedGenes)
      if (input$violin_type_cytokines == 'all'){
        shinyjs::hide('violin_inputCytokines')
        req(data$violin_plot_all)
        data$violin_plot_all
      }
      else{
        shinyjs::show('violin_inputCytokines')
        req(data$violin_plot_input)
        data$violin_plot_input
      }
    }
    else if (input$tabs == 'D'){
      updateSelectizeInput(session, 'boxInputGene', choices = geneList, server = TRUE, selected = box_data$selectedGenes)
      if (input$box_type_cytokines == 'all'){
        shinyjs::hide('box_inputCytokines')
        req(data$box_plot_all)
        data$box_plot_all
      }
      else{
        shinyjs::show('box_inputCytokines')
        req(data$box_plot_input)
        data$box_plot_input
      }
    } else {
      updateSelectizeInput(session, 'umap_gene_input', choices = c(" ", geneList), server = TRUE, selected = umap_data$selectedGenes)
      req(data$plot_res)
      data$plot_res
      # req(data$plot_res2)
      # data$plot_res2
    }
  }
})
