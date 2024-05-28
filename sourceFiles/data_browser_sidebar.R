output$pageStub <- renderUI(fluidPage(
  useShinyjs(),
  tags$head(tags$style(css)),
  titlePanel("Data Browser"),
  sidebarLayout(
    sidebarPanel(
      div(style='overflow-x: auto;',
          radioGroupButtons(inputId = "dataType",
                            choices = c("Gene<br>Expression", "Differential<br>Gene Expression"), 
                            justified = TRUE)),
      tabsetPanel(id = "tabs",
                  tabPanel("Bar Graph", value = "A", 
                           h1("Bar Graph Input"),
                           h1(" "),
                           selectInput("featureInput1", "Choose Cell Type", cellList),
                           selectInput("featureInput2", "Choose Cytokine", newCytokineList),
                  ),
                  tabPanel("Heat Map", value = "B",
                           h1("Heat map Input"),
                           h1(" "),
                           div(style='overflow-x: auto;',
                               radioGroupButtons(inputId = "inputMethod", label = "Method", choices = c("Cell Type", "Cytokine"), justified = TRUE)),
                           selectizeInput("inputGene", "Input Genes",
                                          choices = NULL,
                                          multiple = TRUE),
                           actionButton("sample", "Sample Genes", class = "btn-link"),
                           actionButton("clear", "Clear Genes", class = "btn-link"),
                           HTML('<br><br>'),
                           selectInput("cellType1", "Choose Cell Type", cellList),
                           radioButtons("type_cytokines", "Choose Cytokine Input",
                                        c("All cytokines" = "all",
                                          "Input cytokines" = "inp")),
                           selectInput("inputCytokines", "Choose Cytokine Input", 
                                       choices = newCytokineList,
                                       multiple = TRUE),
                           
                           selectInput("cytokine1", "Choose Cytokine", newCytokineList),
                           radioButtons("type_cellType", "Choose Cell Type Input",
                                        c("All cell types" = "all",
                                          "Input cell types" = "inp")),
                           selectInput("inputCellTypes", "Input Cell Types", 
                                       choices = cellList,
                                       multiple = TRUE),
                           actionButton("submit_celltype", "Submit", class = "btn-block"),
                           actionButton("submit_cytokines", "Submit", class = "btn-block")
                  ),
                  tabPanel("Violin Plot", value = "C",
                           h1("Violin Plot Input"),
                           h1(" "),
                           
                           selectizeInput("violinInputGene", "Input Genes",
                                          choices = NULL,
                                          multiple = TRUE),
                           actionButton("violinSample", "Sample Genes", class = "btn-link"),
                           actionButton("violinClear", "Clear Genes", class = "btn-link"),
                           selectInput("cellType2", "Choose Cell Type", cellList),
                           radioButtons("violin_type_cytokines", "Choose Cytokine Input",
                                        c("All cytokines" = "all",
                                          "Input cytokines" = "inp")),
                           selectInput("violin_inputCytokines", "Choose Cytokine Input", 
                                       choices = cytokineList,
                                       multiple = TRUE),
                           
                           actionButton("violinSubmit", "Submit", class = "btn-block")
                  ),
                  tabPanel("Box Plot", value = "D",
                           h1("Box Plot Input"),
                           h1(" "),
                           selectizeInput("boxInputGene", "Input Genes",
                                          choices = NULL,
                                          multiple = TRUE),
                           actionButton("boxSample", "Sample Genes", class = "btn-link"),
                           actionButton("boxClear", "Clear Genes", class = "btn-link"),

                           selectInput("cellType3", "Choose Cell Type", cellList),

                           radioButtons("box_type_cytokines", "Choose Cytokine Input",
                                        c("All cytokines" = "all",
                                          "Input cytokines" = "inp")),
                           selectInput("box_inputCytokines", "Choose Cytokine Input", 
                                       choices = cytokineList,
                                       multiple = TRUE),

                           actionButton("boxSubmit", "Submit", class = "btn-block")
                  ),
                  tabPanel("UMAP", value = "E",
                            h1("UMAP Input"),
                            h1(" "),
                            selectInput("umap_featureInput1", "Choose Cell Type", c(cellList)),
                            selectInput("umap_featureInput2", "Choose Cytokine", cytokineList),
                            HTML('<hr>'),
                            actionButton("umap_dropdown_btn", label = ' Optional Parameters', icon = icon("chevron-down"),
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
                            hidden(selectizeInput("umap_gene_input", "Gene to Display Expression Level",
                                                  choices = NULL)),
                            # hidden(textInput("umap_gene_input", "Gene to Display Expression Level")),
                            HTML("<br>"),
                            actionButton("umap_submit", "Submit", class = "btn-block")
                  )
      ),
      width = 3
    ),
    mainPanel(
      uiOutput('shinyPlot'),
      HTML('<br><br>'),
      uiOutput('download')
    )
    
  )))

output$shinyPlot <- renderUI(
  tagList(
      div(class = 'barid', style=
              'height: 80vh; 
              width: 110%; 
              overflow-y: auto; 
              position: relative',
              shinycssloaders::withSpinner(plotlyOutput("plot"), type = 5)))
)

output$download <- renderUI({
  if (is.null(input$tabs)) {
  } 
  else if (input$tabs == 'C'){
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