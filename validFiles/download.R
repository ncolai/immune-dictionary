library(shinyjs)
library(shinyWidgets)
library(tidyverse)

# Define the paths
folder_path <- "ref_data/"
extra_file_path <- "dataFiles/lig_seurat.Rda"

# List of data files
data_files <- list.files(folder_path)
full_data_files <- list.files(folder_path)
css <- ".shiny-notification {
          position:fixed;
          top: calc(22%);
          left: calc(25%);
          width: 250px;
}

.item {
  display: flex;
  flex-direction: column;
  align-items: center;
  text-align: center;
}
"

bar <- function(placeholder) {
  output[[placeholder]] <- downloadHandler(
    filename = function(){placeholder},
    content = function(fname){
      foo <- read_rds(paste0(folder_path, placeholder, collapse=""))
      write_rds(foo, fname, compress = "gz")
    }
  )
}

bar_extra <- function() {
  output[["lig_seurat"]] <- downloadHandler(
    filename = function() { "lig_seurat.Rda" },
    content = function(fname) {
      file.copy(extra_file_path, fname)
    }
  )
}

buildDownloadLink <- function(inputData) {
  inputs <- lapply(inputData, function(args) {
    do.call(downloadLink, args)
  })
  for (i in 1:length(full_data_files)) {
    bar(full_data_files[i])
  }
  
  bar_extra()
  
  inputs2 <- list()
  for (i in 1:(3 * length(inputs)-1)) {
    if (i %% 3 == 1) {
      inputs2[[i]] <- img(src = "downloadLinkImg.png", Rd = FALSE, alt = "image", height = 11.5)
    } else if (i %% 3 == 2) {
      inputs2[[i]] <- inputs[[(i + 1) %/% 3]]
    } else {
      inputs2[[i]] <- HTML("<br>")
    }
  }
  
  do.call(div, inputs2)
}

# Create download link list
downloadLinkList <- list()

for (i in 1:length(full_data_files)) {
  downloadLinkList[[i]] <- list(outputId = paste0(full_data_files[[i]]), 
                                label = gsub('ref_data_', '',full_data_files[[i]]))
}

downloadLinkList[[length(full_data_files) + 1]] <- list(
  outputId = "lig_seurat", 
  label = "lig_seurat.Rda"
)

output$record <- renderUI(
  div(id = "downloads", style='display: flex; flex-direction: row; justify-content: space-between;',
      buildDownloadLink(downloadLinkList))
)

output$pageStub <- renderUI(fluidPage(
  useShinyjs(),
  tags$head(tags$style(css)),
  titlePanel("Download"),
  uiOutput("record"),
  div(id = "note",
      style = "margin-top: 20px; padding: 10px; border: 1px solid #ddd; background-color: #f9f9f9;",
      h4("Note: How to Load RDS Files in R"),
      p("How to load RDS files in R, don't worry we can help you! you can use the following code:"),
      code('lig_seurat <- readRDS(paste0("path_to_your_file/"),celltype,".rds",sep = "")'),
      p("Replace ", code('"path_to_your_file.rds"'), " with the actual path to your RDS file."),
      h4("Note: How to Load Multiple RDS Files in R in one Seurat Object"),
      code('lig_seurat <- vector(mode="list", length=15)'),br(),
      code('lig_seurat[[celltype]] <- readRDS(paste0("path_to_your_file/",celltype, ".rds", sep = ""))'),
      p("Replace ", code('"path_to_your_file.rds"'), " with the actual path to your RDS file. Replace ",code('"celltype"'),"with the celltype of interest"),
      h4("lig_seurat.Rda File"),
      p("The file ", code('lig_seurat.Rda'), " contains all cell types combined. This can be useful for comprehensive analysis across multiple cell types."),
      h4("Additional code to see 'celltype', 'sample', 'rep', 'variable', 'value' columns"),
      code('raw_data = lig_seurat[[input$cellType]]'),br(),
      code('plot_df = cbind(t(as.matrix(raw_data@assays[["RNA"]]@data[input$InputGene,])),
                    raw_data@meta.data)'),br(),
      code('plot_df_melt = melt(plot_df, id.vars = names(raw_data@meta.data))'),
      p("Replace ", code('"input$cellType"'), " with the celltype you are interested in. Replace ",code('"input$InputGene"'),"with the gene of interest")
  )
))