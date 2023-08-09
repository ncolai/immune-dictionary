options(shiny.maxRequestSize = 14000 * 1024^5)
options (future.globals.maxSize = 8000 * 1024^5)
library(shiny)
library(dplyr)
library(Seurat)
library(patchwork)
library(googlesheets4)

# put a message in console or server log; note this happens only when the app is started!
cat("uiStub application started...\n")

#drop_load_RData('Web_Signatures/rdata/200527-ligand-allsigs-pcut.RData')
#drop_load_RData('Web_Signatures/rdata/refdata.RData')
#drop_download('Web_Signatures/rdata/refdata.RData', local_path = './refdata.RData', overwrite = TRUE, dtoken = token)
if (!exists('allsigs')){
  cat('Loading data\n')
  load("dataFiles/200527-ligand-allsigs-pcut.RData")
  cat('Loaded data')
}

# table <- "responses"
# 
# saveData <- function(data) {
#   # The data must be a dataframe rather than a named vector
#   data <- data %>% as.list() %>% data.frame()
#   # Add the data as a new row
#   sheet_append(SHEET_ID, data)
# }

css <- ".shiny-notification {
          position:fixed;
          top: calc(20%);
          left: calc(10%);
        }"

loginpage <- 
  tagList(             # a single-output stub ui basically lets you
    fluidPage(                                  #     move the ui into the server function
      fluidRow(
        tags$head(tags$style(css)),
        column(12,
               includeHTML("htmlFiles/logIn.html")
        )
      )
    ),
    div(id = "login_panel",style = "width: 90%; padding: 20px;",
     wellPanel(textInput("userName", "Please enter your email here to access the website: "),
     br(),actionButton("login", "Log in")))
  )
  
#   tagList(div(id = "login", style = "width: 80%; padding: 20px;",
#          wellPanel(textInput("userName", "Username"),
#          passwordInput("passwd", "Password"),
#          br(),actionButton("login", "Log in")),
#          
# ))
  
#   div(id = "loginpage", style = "width: 500px; max-width: 100%; margin: 0 auto; padding: 20px;",
#      wellPanel(
#        tags$h2("LOG IN", class = "text-center", style = "padding-top: 0;color:#333; font-weight:600;"),
#        textInput("userName", placeholder="Username", label = tagList(icon("user"), "Username")),
#        passwordInput("passwd", placeholder="Password", label = tagList(icon("unlock-alt"), "Password")),
#        br(),
#        div(
#          style = "text-align: center;",
#          actionButton("login", "SIGN IN", style = "color: white; background-color:#3c8dbc;
#                      padding: 10px 15px; width: 150px; cursor: pointer;
#                      font-size: 18px; font-weight: 600;"),
#        ))
# )

# credentials = data.frame(
#   username_id = c("myuser", "myuser1"),
#   passod   = sapply(c("mypass", "mypass1"),password_store),
#   permission  = c("basic", "advanced"), 
#   stringsAsFactors = F
# )

login = FALSE
first = TRUE
USER <- reactiveValues(login = login, first = first)

ui <- uiOutput("uiStub")                                # single-output stub ui

isValidEmail <- function(x) {
  grepl("\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>", as.character(x), 
        ignore.case=TRUE)
}

server <- function(input, output, session) {
  cat("Session started.\n")                               # this prints when a session starts
  onSessionEnded(function() {cat("Session ended.\n\n")})  # this prints when a session ends
  
  
  
  # observeEvent(input$login,{
  #   if (USER$login == FALSE) {
  #       if (input$userName != '' && isValidEmail(input$userName)) {
  #         USER$login <- TRUE
  #         name <- input$userName %>% as.list()
  #         time <- c(as.character(Sys.time()))
  #         zone <- c(as.character(Sys.timezone()))
  #         data <- data.frame(name, time, zone)
  #         
  #         # drive_auth(path = "client_secret.json")
  #         gs4_auth(path = "client_secret.json")
  #         
  #         sheet_append(ss = 'https://docs.google.com/spreadsheets/d/11yyHKU8VkvKeK8M79iEYk_joqtEV5hBDU_s1X1-UT3c/edit?usp=sharing',
  #                      data = data,
  #                      sheet = 1)
  #       }
  #       else{
  #         showNotification('Please enter a valid email.', type = 'error')
  #       }
  #   }   
  # }
               
  # )
  
  # tags$script('
  #         $(document).ready(function(){
  #         var d = new Date();
  #         var target = $("#clientTime");
  #         target.val(d.toLocaleString());
  #         target.trigger("change");
  #         });
  #         ')
  # textInput("clientTime", "Client Time", value = "")
  
  # build menu; same on all pages
  output$uiStub <- renderUI({
    if (USER$login == FALSE){
      tagList(             # a single-output stub ui basically lets you
      fluidPage(                                  #     move the ui into the server function
        fluidRow(
          column(12,
                 includeHTML("htmlFiles/mainPage.html")
          )
        ),
        uiOutput("pageStub")                     # loaded server code should render the
        
        )                                           #    rest of the page to this output$
      )
    }
    else {
      loginpage
    }})
  
  # load server code for page specified in URL
  validFiles = c("home.R",                             # valid files must be hardcoded here
                 "paper.R",
                 "data_browser.R",
                 "irea.R",
                 "download.R",
                 "instructions.R",
                 "contact.R")                     #    for security (use all lower-case
                                                  #    names to prevent Unix case problems)
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
  source(file.path("validFiles", fname), local=TRUE)                            # load and run server code for this page
}
# Run the application
shinyApp(ui = ui, server = server)


