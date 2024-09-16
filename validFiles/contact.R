library(mailR)
library(rJava)

css <- ".shiny-notification {
          position:fixed;
          top: calc(30%);
          left: calc(35%);
        }"

output$pageStub <- renderUI(fluidPage(
  fluidRow(
    tags$head(tags$style(css)),
    HTML('<br>'),
    p('If you have any comments or concerns, please do not hesitate to reach out to us at: ImmuneDictionary@gmail.com')
)))
