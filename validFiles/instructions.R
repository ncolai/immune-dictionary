library(shinyjs)
library(shinyWidgets)
library(rmarkdown)

css <- "
.shiny-notification {
          position:fixed;
          top: calc(22%);
          left: calc(25%);
          width: 250px;
        }
"

output$pageStub <- renderUI(fluidPage(
  useShinyjs(),
  tags$head(tags$style(css)),
  titlePanel("Instructions"),
  sidebarLayout(
    sidebarPanel(
      includeHTML(render('www/tableOfContent.Rmd'))
    ),
    mainPanel (
      uiOutput("header")
    )
  )
))

output$header <- renderUI(tagList(
  includeHTML(render('www/instructionContent.Rmd'))
))