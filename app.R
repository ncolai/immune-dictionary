options(shiny.maxRequestSize = 7000 * 1024^5)
options (future.globals.maxSize = 4000 * 1024^5)
library(shiny)
library(dplyr)
library(patchwork)
library(googlesheets4)
library(fst)
library(qs)
library(Seurat)

# put a message in console or server log; note this happens only when the app is started!
cat("uiStub application started...\n")


css <- ".shiny-notification {
          position:fixed;
          top: calc(20%);
          left: calc(10%);
        }"

ui <- fluidPage(
  tags$head(includeHTML("htmlFiles/google_analytics.html")),
  uiOutput("uiStub")                               # single-output stub ui
)

server <- function(input, output, session) {
  cat("Session started.\n")                               # this prints when a session starts
  onSessionEnded(function() {cat("Session ended.\n\n")})  # this prints when a session ends
  
  # build menu; same on all pages
  output$uiStub <- renderUI({
    tagList(             # a single-output stub ui basically lets you
      fluidPage(            # move the ui into the server function
        includeCSS("htmlFiles/styles.css"),
        includeHTML("htmlFiles/mainPage.html"),
        # fluidRow(
        #   column(12, includeHTML("htmlFiles/mainPage.html"))
        # ),
        uiOutput("pageStub")                  # loaded server code should render the
      )                                           #    rest of the page to this output$
    )
  })
  
  # load server code for page specified in URL
  validFiles = c("home.R",
                 "data_browser.R",
                 "irea.R",
                 "ref_card.R",
                 "download.R",
                 "instructions.R",
                 "contact.R")                     #    for security (use all lower-case
  fname = isolate(session$clientData$url_search)       # isolate() deals with reactive context
  if(nchar(fname)==0) { fname = "?home" }              # blank means home page
  fname = paste0(substr(fname, 2, nchar(fname)), ".R") # remove leading "?", add ".R"
  
  
  if(!fname %in% validFiles){                          # is that one of our files?
    output$pageStub <- renderUI(tagList(              # 404 if no file with that name
      fluidRow(
        column(5,
               HTML("<h2>404 Not Found Error:</h2><p>That URL doesn't exist. Use the",
                    "menu above to navigate to the page you were looking for.</p>")
        )
      )
    ))
    return()    # to prevent a "file not found" error on the next line after a 404 error
  }
  source(file.path("validFiles", fname), local=TRUE)    # load and run server code for this page
}
shinyApp(ui = ui, server = server)
