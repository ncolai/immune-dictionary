library(mailR)

css <- ".shiny-notification {
          position:fixed;
          top: calc(30%);
          left: calc(35%);
        }"

output$pageStub <- renderUI(fluidPage(
  fluidRow(
    tags$head(tags$style(css)),
    HTML('<br>'),
    div(style=
          'padding: 0;
                        margin: 0;
                        color: black;
                        font-weight: bold;
                        font-size: 1em;
                        white-space: normal;
                        width: 100%;
                        text-align: left;
                        margin-bottom: 5px;',
        HTML('Contact Form')
    ),
    div(id = "login",
        wellPanel(title = "Contact Us", 
                  # textInput("to", label = "To:", placeholder = "To:"),
                  textInput("name", label = "Name: *", placeholder = "Name:"),
                  textInput("from", label = "Email: *", placeholder = "Email:"),
                  textInput("sub","Subject: *", placeholder = 'Subject:'),
                  textAreaInput(inputId="message", label = "Message: *", rows = 15, placeholder = 'Message here...', resize = "vertical"),
                  actionButton("mailButton",label = "Send mail") 
        )
    ))
))

isValidEmail <- function(x) {
  grepl("\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>", as.character(x), 
        ignore.case=TRUE)
}


observeEvent(input$mailButton,{
  if (input$from == '' | input$sub == '' | input$message == '' | input$name == ''){
    showNotification('Please fill out all required fields.', type = 'error')
  }
  else{
    if (isValidEmail(input$from)){
      print('sending mail')
      send.mail(from = "immunedictionary@gmail.com",
                to = c("byao@broadinstitute.org", "angcui@mit.edu"),
                subject = paste('[Immune Dictionary]', input$sub),
                body = paste('Name: ', input$name, '\n\nEmail: ', input$from, '\n\n', 'Subject: ', input$sub, '\n\nMessage: \n', input$message, sep = ''),
                smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "immunedictionary@gmail.com", passwd = "IREA2021", ssl = TRUE),
                authenticate = TRUE,
                send = TRUE)
      showNotification('Your message has been sent! We will contact you shortly.')
      print('sent email!')
    }
    else{
      showNotification('Please enter a valid email.', type = 'error')
    }
  }
  
  
  # send.mail(from = input$from,
  #             to = input$to,
  #             subject = input$sub,
  #             body = input$msg,
  #             smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "myemail@gmail.com", passwd = "mypasword", ssl = TRUE),
  #             authenticate = TRUE,
  #             html = TRUE,
  #             send = TRUE)
})
