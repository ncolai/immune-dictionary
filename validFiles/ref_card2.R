library(shinyjs)
library(shinyWidgets)

css <- "
.shiny-notification {
          position:fixed;
          top: calc(22%);
          left: calc(25%);
          width: 250px;
        }
"

cellList <- readLines("sourceFiles/lig_seurat_data.txt")
cellList <- cellList[cellList != "eTAC"]
cellList <- cellList[cellList != "ILC"]
cellList <- cellList[cellList != "Mast_cell"]
cellList <- c(cellList, "Langerhans")
#cellList <- sort(cellList)

output$pageStub <- renderUI(tagList(
  useShinyjs(),
  tags$head(tags$style(css)),
  titlePanel("Reference Cards for Immune Cell Polarization States"),
  HTML("Figure 3 and its associated figures of the publication \"Dictionary of immune responses to cytokines at single-cell resolution\". <br><br>"),
  sidebarLayout(
    sidebarPanel(
      selectInput("ref_card_cell", "Cell Type", cellList),
      width = 4
    ), mainPanel(
      uiOutput("pdfOutput")
    )
  )
))

output$pdfOutput <- renderUI(tagList(
  tags$style("display: flex; justify-content: center; align-items: center; height: 100vh;"),
  tags$img(
    src = paste0("ref_cards/", input$ref_card_cell, ".png"),
    style = "width: 820px; height: 1040px;",
    alt=paste0(input$ref_card_cell, " Reference Card")
  )
))