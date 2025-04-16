first_visit_tab <- reactiveValues(A = TRUE,
                                  B = TRUE,
                                  C = TRUE,
                                  D = TRUE,
                                  E = TRUE)

observe({
  req(input$inputGene)
  if (is.null(input$tabs)) {
  } else if (input$tabs == 'B' && first_visit_tab$B){
    #click('submit_celltype')
    #click('submit_cytokines')
    data$plot_all1 <- example_plots$heatmap_plot_1
    data$plot_all2 <- example_plots$heatmap_plot_2
    #heatmap_data$selectedGenes <- c("Il4i1", "St7", "Isg15")
    first_visit_tab$B = FALSE
  }
})

observe({
  req(input$violinInputGene)
  if (is.null(input$tabs)) {
  } else if (input$tabs == 'C' && first_visit_tab$C){
    showNotification('Sample Graph is Loading.', type = 'message')
    data$violin_plot_all <- example_plots$violin_plot_example
    first_visit_tab$C = FALSE
  }
})

observe({
  req(input$boxInputGene)
  if (is.null(input$tabs)) {
  } else if (input$tabs == 'D' && first_visit_tab$D){
    showNotification('Sample Graph is Loading.', type = 'message')
    data$box_plot_all <- example_plots$box_plot_example
    first_visit_tab$D = FALSE
  }
})

observe({
  req(input$umap_featureInput1)
  if (is.null(input$tabs)) {
  } else if (input$tabs == 'E' && first_visit_tab$E){
    showNotification('Sample Graph is Loading.', type = 'message')
    data$plot_res <- example_plots$umap_plot_example
    first_visit_tab$E = FALSE
  }
})