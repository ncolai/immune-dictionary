output$pageStub <- renderUI(fluidPage(
  useShinyjs(),
  tags$head(tags$style(css)),
  titlePanel("Cytokine DEG"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(id = "tabs",
                  tabPanel("Bar Graph", value = "A", 
                           h1("Bar Graph Input"),
                           h1(" "),
                           selectInput("featureInput1", "Choose Cell Type", cellList),
                           selectInput("featureInput2", "Choose Cytokine", cytokineList),
                  ),
                  tabPanel("Heat Map", value = "B",
                           h1("Heat map Input"),
                           h1(" "),
                           div(style='overflow-x: auto;',
                               radioGroupButtons(inputId = "inputMethod", label = "Method", choices = c("Cell Type", "Cytokine"), justified = TRUE)),
                           textAreaInput(inputId="inputGene", label = "Input Genes", rows = 15, placeholder = 'Gene1\nGene2\n...', resize = "vertical"),
                           actionButton("sample", "Sample Genes", class = "btn-link"),
                           actionButton("clear", "Clear Genes", class = "btn-link"),
                           HTML('<br><br>'),
                           selectInput("cellType1", "Choose Cell Type", cellList),
                           radioButtons("type_cytokines", "Choose Cytokine Input",
                                        c("All cytokines" = "all",
                                          "Input cytokines" = "inp")),
                           
                           textAreaInput(inputId="inputCytokines", label = "Input Cytokines", rows = 15, placeholder = 'Cytokine1\nCytokine2\n...', resize = "vertical"),
                           
                           selectInput("cytokine1", "Choose Cytokine", cytokineList),
                           radioButtons("type_cellType", "Choose Cell Type Input",
                                        c("All cell types" = "all",
                                          "Input cell types" = "inp")),
                           textAreaInput(inputId="inputCellTypes", label = "Input Cell Types", rows = 15, placeholder = 'Cell Type1\nCell Type2\n...', resize = "vertical"),
                           actionButton("submit_celltype", "Submit", class = "btn-block"),
                           actionButton("submit_cytokines", "Submit", class = "btn-block")
                  ),
                  tabPanel("Violin Plot", value = "C",
                           h1("Violin Plot Input"),
                           h1(" "),
                           textAreaInput(inputId="violinInputGene", label = "Input Genes", rows = 15, placeholder = 'Gene1\nGene2\n...', resize = "vertical"),
                           actionButton("violinSample", "Sample Genes", class = "btn-link"),
                           actionButton("violinClear", "Clear Genes", class = "btn-link"),
                           
                           selectInput("cellType2", "Choose Cell Type", cellList),
                           
                           radioButtons("violin_type_cytokines", "Choose Cytokine Input",
                                        c("All cytokines" = "all",
                                          "Input cytokines" = "inp")),
                           textAreaInput(inputId="violin_inputCytokines", label = "Input Cytokines", rows = 15, placeholder = 'Cytokine1\nCytokine2\n...', resize = "vertical"),
                           
                           actionButton("violinSubmit", "Submit", class = "btn-block")
                  ),
                  tabPanel("Box Plot", value = "D",
                           h1("Box Plot Input"),
                           h1(" "),
                           textAreaInput(inputId="boxInputGene", label = "Input Genes", rows = 15, placeholder = 'Gene1\nGene2\n...', resize = "vertical"),
                           actionButton("boxSample", "Sample Genes", class = "btn-link"),
                           actionButton("boxClear", "Clear Genes", class = "btn-link"),

                           selectInput("cellType3", "Choose Cell Type", cellList),

                           radioButtons("box_type_cytokines", "Choose Cytokine Input",
                                        c("All cytokines" = "all",
                                          "Input cytokines" = "inp")),
                           textAreaInput(inputId="box_inputCytokines", label = "Input Cytokines", rows = 15, placeholder = 'Cytokine1\nCytokine2\n...', resize = "vertical"),

                           actionButton("boxSubmit", "Submit", class = "btn-block")
                  )
      ),
      width = 3
    ),
    mainPanel(
      div(class = 'barid', style=
            'height: 80vh; 
            width: 115vh; 
            overflow-y: auto; 
            position: relative',
          shinycssloaders::withSpinner(plotlyOutput("plot"), type = 5)),
      HTML('<br><br>'),
      uiOutput('download')
      
      #plotlyOutput("heatmap1", height = 500),
    )
    
  )))

output$download <- renderUI({
  if (input$tabs == 'C'){
    if (input$violin_type_cytokines == 'all'){
      uiOutput("download_violin_all")
    } else {
      uiOutput("download_violin_inp")
    }
  }
  else if (input$tabs == 'D'){
    if (input$box_type_cytokines == 'all'){
      uiOutput("download_box_all")
    } else {
      uiOutput("download_box_inp")
    }
  }
})