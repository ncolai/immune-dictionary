library(plotly)
library(tidyverse)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders)
library(reshape2)

cellList <- c(" ", readLines("sourceFiles/lig_seurat_data.txt"))
cytokineList <- c(" ", readLines("sourceFiles/cytokine_list.txt"), "")
lig_seurat <- vector(mode="list", length=15)
names(lig_seurat) <- readLines("sourceFiles/lig_seurat_data.txt")
source('sourceFiles/cytokine_sidebar.R', local = TRUE)
source('sourceFiles/cytokine_heat_map.R', local = TRUE)
source('sourceFiles/cytokine_violin_plot.R', local = TRUE)
source('sourceFiles/cytokine_box_plot.R', local = TRUE)


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
                       box_input_data = NULL)

"%ni%" <- Negate("%in%")
valid_comp <- function(gene_list, valid_genes) {
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
  if (input$tabs == 'A'){
    if (input$featureInput1 == " ") return(NULL)
    else if (input$featureInput2 == " ") return(NULL)
    else{
      subset_data <- subset(allsigs, celltype == input$featureInput1)
      final_data <- subset(subset_data, sample == input$featureInput2)
      gene_column <- final_data %>% pull(gene)
      log_column <- final_data %>% pull(avg_logFC)
      
      y <- factor(gene_column)
      x <- log_column
      plot_height <- 500 + 15*nrow(final_data)
      
      fig <- plot_ly(
        y = reorder(y,x),
        x = x,
        name = "DEG Cytokine Expression",
        type = "bar",
        text = paste0("p-value: ",final_data %>% pull(p_val),sep=""),
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
  else if (input$tabs == 'B'){
    if (input$inputMethod == "Cell Type"){
      show('cellType1')
      # show('inputCytokines')
      show('submit_celltype')
      show('type_cytokines')
      hide('cytokine1')
      hide('inputCellTypes')
      hide('submit_cytokines')
      hide('type_cellType')
      
      if (input$type_cytokines == 'all'){
        hide('inputCytokines')
        req(data$plot_all1)
        data$plot_all1
      }
      else{
        show('inputCytokines')
        req(data$plot_input1)
        data$plot_input1
      }
    }
    else {
      hide('cellType1')
      hide('inputCytokines')
      hide('submit_celltype')
      hide('type_cytokines')
      show('cytokine1')
      #show('inputCellTypes')
      show('submit_cytokines')
      show('type_cellType')
      
      if (input$type_cellType == 'all'){
        hide('inputCellTypes')
        req(data$plot_all2)
        data$plot_all2
      }
      else{
        show('inputCellTypes')
        req(data$plot_input2)
        data$plot_input2
      }
    }
  }
  else if (input$tabs == 'C'){
    if (input$violin_type_cytokines == 'all'){
      hide('violin_inputCytokines')
      req(data$violin_plot_all)
      data$violin_plot_all
    }
    else{
      show('violin_inputCytokines')
      req(data$violin_plot_input)
      data$violin_plot_input
    }
  }
  else if (input$tabs == 'D'){
    if (input$box_type_cytokines == 'all'){
      hide('box_inputCytokines')
      req(data$box_plot_all)
      data$box_plot_all
    }
    else{
      show('box_inputCytokines')
      req(data$box_plot_input)
      data$box_plot_input
    }
  }
}
)
