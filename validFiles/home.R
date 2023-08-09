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
      
      .paper-link {
        position:relative;
        top: 200px;
      }"),

  div(
    id = "title",
    "Immune Dictionary"
  ),
  div(
    id = "description",
    "Cytokines mediate cell-cell communication in the immune system and represent important therapeutic targets. While there have been in-depth studies of individual cytokines, we lack a global view of the responses of each immune cell type to each cytokine. To address this gap we created the Immune Dictionary – a compendium of single-cell transcriptomic profiles of over 20 cell types in response to each of 86 cytokines in murine lymph nodes in vivo. Based on the dictionary, we developed companion software, Immune Response Enrichment Analysis (IREA), for assessing cytokine activities and immune cell polarization from transcriptomic data, and applied it to reveal cytokine networks in tumors following immune checkpoint blockade therapy. Our dictionary generates new hypotheses for cytokine functions, illuminates pleiotropic effects of cytokines, expands our knowledge of activation states in each immune cell type, and provides a framework to deduce the roles of specific cytokines and cell-cell communication networks in any immune response."
  ),
  p( style="text-align:center",
    a(
      "This will be the link to the paper",
      target = "_blank",
      href = "https://www.google.com/",
      class = "paper-link"
    )
  )
  # Maybe transfer into another file when UI is figured out
  
  # includeHTML("htmlFiles/home.html")
))
