library(Seurat)
library(ggplot2)
library(patchwork)
library(plotly)
library(shinyjs)
library(cowplot)

cellList <- readLines("sourceFiles/lig_seurat_data.txt")
lig_seurat_cells <- vector(mode="list", length=15)
names(lig_seurat_cells) <- cellList
cytokineList <- c(" ", readLines("sourceFiles/cytokine_list.txt"), "")

output$pageStub <- renderUI(fluidPage(
  useShinyjs(),
  titlePanel("UMAP"),
  sidebarLayout(
    sidebarPanel(
    fluidRow(
      selectInput("featureInput1", "Choose Cell Type", c(" ", cellList)),
      selectInput("featureInput2", "Choose Cytokine", cytokineList),
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
      hidden(textInput("gene_input", "Gene to Display Expression Level")),
      HTML("<br>"),
      actionButton("submit", "Submit", class = "btn-block"),
      )
    )
  ,
  mainPanel(
    shinycssloaders::withSpinner(plotlyOutput("plot", height = 600, width = 600), type = 5),
    shinycssloaders::withSpinner(plotlyOutput("plot2", height = 600, width = 600), type = 5)
  )
  )
  
))

valid_comp <- function(gene, valid_genes) {
  if (gene %in% valid_genes){
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}

data_vals <- reactiveValues(plot_res1 = NULL, plot_res2 = NULL)

observeEvent(input$dropdown_btn,{
  # hide/show optional parameters on click
  toggle('gene_input')
})

observeEvent(input$submit, {
  if (input$featureInput1 != ' ' & input$featureInput2 != ' '){
    input_celltype = input$featureInput1
    if (is.null(lig_seurat_cells[[input_celltype]])){
      print('loading data')
      file_name = paste0('dataFiles/celltype/200417-ligands-alldata-seurat-p3-', input$featureInput1, '.RDS', sep = '')
      lig_seurat_cells[[input_celltype]] <<- readRDS(file_name)
    }
    
    if (input$featureInput2 == 'ALL'){
      baseplot <- DimPlot(lig_seurat_cells[[input_celltype]], reduction = "umap", pt.size = 1, group.by = "sample", label = FALSE)
    }
    else{
      names <- readLines("sourceFiles/cytokine_list.txt")
      values <- rep('red', times = length(names))
      inp_cytokine <- list()

      for (i in seq_along(names)) {
        inp_cytokine[[names[[i]]]] <- setNames(values[i], names[i])
      }

      baseplot <- DimPlot(lig_seurat_cells[[input_celltype]], reduction = "umap", pt.size = 1, group.by = "sample", cols = inp_cytokine[[input$featureInput2]], order = input$featureInput2, label = FALSE)
    }
    data_vals$plot_res1 <- baseplot + ggtitle(paste0(input$featureInput1, " UMAP", sep = '')) + NoLegend() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                              axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                                                                              axis.title.x=element_blank(),
                                                                                                              axis.title.y=element_blank(),legend.position="none",
                                                                                                              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                                                                                              panel.grid.minor=element_blank(),plot.background=element_blank())

    if (input$gene_input != '' & valid_comp(input$gene_input, allsigs$gene)){
      # check if gene is valid
      data_vals$plot_res2 <- FeaturePlot(lig_seurat_cells[[input_celltype]], features = c(input$gene_input)) + NoLegend() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                                             axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                                                                                             axis.title.x=element_blank(),
                                                                                                                             axis.title.y=element_blank(),legend.position="none",
                                                                                                                             panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                                                                                                             panel.grid.minor=element_blank(),plot.background=element_blank())
    }
    # plot2 <- FeaturePlot(lig_seurat_cells[[input_celltype]], features = c("Ifitm3"))
    # (baseplot + NoLegend()) + plot2 + plot_layout(ncol=1)
    # (baseplot + ggtitle(paste0(input$featureInput1, " UMAP", sep = '')) + NoLegend() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
    #                               axis.text.y=element_blank(),axis.ticks=element_blank(),
    #                               axis.title.x=element_blank(),
    #                               axis.title.y=element_blank(),legend.position="none",
    #                               panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
    #                               panel.grid.minor=element_blank(),plot.background=element_blank())) + plot2

  }
})

# observeEvent(input$submit, {
#   if (input$featureInput1 != ' ' & input$featureInput2 != ' '){
#     input_celltype = input$featureInput1
#     if (is.null(lig_seurat_cells[[input_celltype]])){
#       print('loading data')
#       file_name = paste0('celltype/200417-ligands-alldata-seurat-p3-', input$featureInput1, '.RDS', sep = '')
#       lig_seurat_cells[[input_celltype]] <<- readRDS(file_name)
#     }
#     
#     object_lig = lig_seurat_cells[[input_celltype]]
#     
#     # reduction <- "umap" %||% DefaultDimReduc(object = object)
#     cells <- colnames(x = object_lig)
#     dims <- c(1,2)
#     data <- Embeddings(object = object_lig[["umap"]])[cells, dims]
#     print('hi')
#     data <- as.data.frame(x = data)
#     dims <- paste0(Key(object = object_lig[["umap"]]), dims)
#     # dims <- paste0(Key(object = object_lig[["umap"]]), dims)
#     # print('hello')
#     object_lig[['ident']] <- Idents(object = object_lig)
#     group.by <- "sample"
#     data <- cbind(data, object_lig[[group.by]][cells, , drop = FALSE])
#     group.by <- colnames(x = data)[3:ncol(x = data)]
#     # print('gud')
#     for (group in group.by) {
#       print(group)
#       if (!is.factor(x = data[, group])) {
#         data[, group] <- factor(x = data[, group])
#       }
#     }
#     
#     col.by <- "sample"
#     pt.size <- 1
#     order <- "41BBL"
#     order <- rev(x = c(
#       order,
#       setdiff(x = unique(x = data[, col.by]), y = order)
#     ))
#     data[, col.by] <- factor(x = data[, col.by], levels = order)
#     new.order <- order(x = data[, col.by])
#     data <- data[new.order, ]
#     if (length(x = pt.size) == length(x = new.order)) {
#       pt.size <- pt.size[new.order]
#     }
#     
#     plot <- ggplot(data = data)
#     plot <- plot + geom_point(
#       mapping = aes_string(
#         x = dims[1],
#         y = dims[2],
#         col = paste0("`", "sample", "`"),
#         shape = NULL,
#         alpha = NULL
#       ),
#       size = 1
#     )
#     
#     plot <- plot +
#       guides(color = guide_legend(override.aes = list(size = 3))) +
#       labs(color = NULL, title = 'UMAP') +
#       CenterTitle()
#     
#     cols <- c('41BBL' = 'red')
#     # colors <- DiscretePalette(length(unique(data[["sample"]])), palette = c('41BBL' = 'red'))
#     # colors <- DiscretePalette(length(unique(data[["sample"]])), palette = cols)
#     # scale <- scale_color_manual(values = colors, na.value = "grey50")
#     # (values = cols, na.value = 'grey50')
#     scale <- scale_color_manual(values = cols, na.value = "grey50")
#     # 
#     # plot <- plot + theme_cowplot()
#     
#     data_vals$plot_res1 <- plot + scale + theme_cowplot() + NoLegend()
#     # 
#     # object[['ident']] <- Idents(object = object)
#     # orig.groups <- "sample"
#     # group.by <- "sample" %||% 'ident'
#     # data <- cbind(data, object[[group.by]][cells, , drop = FALSE])
#     # group.by <- colnames(x = data)[3:ncol(x = data)]
#     # for (group in group.by) {
#     #   if (!is.factor(x = data[, group])) {
#     #     data[, group] <- factor(x = data[, group])
#     #   }
#     # }
#     
#   }
# })

output$plot <- renderPlotly({
  # req(input$featureInput1)
  # req(input$featureInput2)
  data_vals$plot_res1
})

output$plot2 <- renderPlotly({
  # req(input$featureInput1)
  # req(input$featureInput2)
  # req(input$gene_input)
  data_vals$plot_res2
})

