library(shinyjs)
library(shinyWidgets)
library(tidyverse)

# get files from the folder
folder_path <- "ref_data/"
data_files <- list.files(folder_path)
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

# email <- reactiveValues(login = FALSE)

bar <- function(placeholder) {
  output[[placeholder]] <- downloadHandler(
    filename = function(){placeholder},
    content = function(fname){
      foo <- read_rds(paste0(folder_path, placeholder, collapse=""))
      write_rds(foo, fname, compress = "gz")
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

folder_path <- "ref_data/"
full_data_files <- list.files(folder_path)

downloadLinkList <- list()

for (i in 1:length(full_data_files)) {
  downloadLinkList[[i]] <- list(outputId = paste0(full_data_files[[i]]), 
                                label = gsub('ref_data_', '',full_data_files[[i]]))
}

# isValidEmail <- function(x) {
#   grepl("\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>", as.character(x),
#         ignore.case=TRUE)
# }

output$record <- renderUI(
  # if (email$login) {
  div(id = "downloads", style='display: flex; flex-direction: row; justify-content: space-between;',
      buildDownloadLink(downloadLinkList))
  # }
  # else {
  #   div(id = "login_panel",style = "width: 90%; padding: 20px;",
  #       wellPanel(textInput("userName", "Please enter your email here to access the website: "),
  #                 br(),actionButton("loginButton", "Log in")))
  # }
)

# observeEvent(input$loginButton,{
#   if (isValidEmail(input$userName)) {
#     email$login = TRUE
#     name <- input$userName %>% as.list()
#     time <- c(as.character(Sys.time()))
#     zone <- c(as.character(Sys.timezone()))
#     data <- data.frame(name, time, zone)
# 
#     # drive_auth(path = "client_secret.json")
#     gs4_auth(path = "client_secret.json")
# 
#     sheet_append(ss = 'https://docs.google.com/spreadsheets/d/11yyHKU8VkvKeK8M79iEYk_joqtEV5hBDU_s1X1-UT3c/edit?usp=sharing',
#                  data = data,
#                  sheet = 1)
#   } else {
#     showNotification("Please enter a valid email address", type = "error")
#   }
# })

output$pageStub <- renderUI(fluidPage(
  useShinyjs(),
  tags$head(tags$style(css)),
  titlePanel("Download"),
  uiOutput("record")
))
