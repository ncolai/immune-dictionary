output$pageStub <- renderUI(tagList(
  tags$style("
      #title {
        text-align: center;
        position: relative;
        top: 0px;
        font-size: 32px;
      }

      #description {
        text-align: center;
        position: relative;
        top: 50px;
        font-size: 16px;
      }
      
      .center {
        font-size: 8px;
        display: block;
        margin-top: 3%;
        margin-left: auto;
        margin-right: auto;
        width: 75%;
        height: 50%;
      }"),
  HTML('<br>'),
  div(
    id = "title",
    strong("Welcome to the Immune Dictionary Portal!")
  ),
  div(
    id = "description",
    "The Immune Dictionary is a compendium of single-cell transcriptomic profiles of over 17 cell types in response to each of 86 cytokines."
  ),
  HTML('<br><br><br>'),
  p( style="text-align:center",
     "More Details on the Publication:",
     HTML("<br>Cui, A. et al., Dictionary of immune responses to cytokines at single-cell resolution. <b>Nature</b>. "),
     a(
       "https://doi.org/10.1038/s41586-023-06816-9",
       target = "_blank",
       href = "https://doi.org/10.1038/s41586-023-06816-9",
       class = "paper-link"
     ),
     "(2023)"
  ),
  HTML(''),
  img(
    class = 'center',
    src = "coverImage.png", 
    alt = "Cytokines Cover Image",
  ),
  div(
    id = "description",
    "Immune Response Enrichment Analysis (IREA) is the companion software for the Immune Dictionary. IREA assesses cytokine activities, immune cell polarization, and cell-cell communication from user transcriptomic data."
  ),
  HTML('<br><br>'),
  img(
    class = 'center',
    src = "IREACoverImage.png", 
    alt = "IREA Cover Image",
  ),
))
