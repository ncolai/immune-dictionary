observeEvent(input$violinSample,{
  updateSelectizeInput(session,
                    "violinInputGene",
                    choices=geneList,
                    server=TRUE,
                    selected = c("Il4i1", "St7", "Isg15")
  )
})

observeEvent(input$violinClear,{
  updateSelectInput(session, "violinInputGene", selected = "")
})

output$download_violin_all <- renderUI({
  req(data$violin_plot_all)
  tagList(downloadButton('download_violins_all_csv', 'Download Violin Plot as CSV'),
          downloadButton('download_violins_all_pdf', 'Download Violin Plot as PDF'),
          HTML("<br><br><br>"))
})

output$download_violin_inp <- renderUI({
  req(data$violin_plot_input)
  tagList(downloadButton('download_violins_inp_csv', 'Download Violin Plot as CSV'),
          downloadButton('download_violins_inp_pdf', 'Download Violin Plot as PDF'),
          HTML("<br><br><br>"))
})

output$download_violins_all_csv <- downloadHandler(
  filename = function(){"violin_plot_all.csv"},
  content = function(fname){
    write.csv(data$violin_all_data, fname)
  }
)

output$download_violins_inp_csv <- downloadHandler(
  filename = function(){"violin_plot_inp.csv"},
  content = function(fname){
    write.csv(data$violin_input_data, fname)
  }
)

# final_plot_violin <- function(num) {
#   if (num == 1) {
#     req(data$violin_plot_all)
#     plot_all <- data$violin_plot_all
#     plot_all
#   } else if (num == 2) {
#     req(data$violin_plot_inp)
#     plot_inp <- data$violin_plot_inp
#     plot_inp
#   }
# }

output$download_violins_all_pdf <- downloadHandler(
  filename = function(){paste("violin_plot_all ",Sys.Date(),".pdf",sep="")},
  content = function(file){
    saveWidget(data$violin_plot_all, 'temp.html', selfcontained = FALSE)
    webshot2::webshot('temp.html', file = file, zoom = 5)
    unlink('temp.html')
  }
)

output$download_violins_inp_pdf <- downloadHandler(
  filename = function(){paste("violin_plot_inp ",Sys.Date(),".pdf",sep="")},
  content = function(file){
    saveWidget(data$violin_plot_inp, 'temp.html', selfcontained = FALSE)
    webshot2::webshot('temp.html', file = file, zoom = 5)
    unlink('temp.html')
  }
)

# output$download_violins_all_pdf <- downloadHandler(
#   filename = "violin_plot_all.pdf",
#   content = function(file) {
#     plot_to_save <- final_plot_violin(1)
#     save_image(plot_to_save, file = file, scale = 5)
#   }) 
# 
# output$download_violins_inp_pdf <- downloadHandler(
#   filename = "violin_plot_inp.pdf",
#   content = function(file) {
#     plot_to_save <- final_plot_violin(2)
#     save_image(plot_to_save, file = file, scale = 5)
#   }) 

submit_btn_celltype2 <- observeEvent(input$violinSubmit,{
  suppressWarnings({
  start_time <- Sys.time()
  if (input$cellType2 == " ") {
    showNotification('Please select a cell type.', type = 'error')
  }
  else if (length(input$violinInputGene) == 0) {
    showNotification('Please enter in genes.', type = 'error')
  }
  else if (input$violin_type_cytokines == 'inp' & length(input$violin_inputCytokines) > 0){
    var_name <- paste0("ft2_", input$cellType2)
    if (exists(var_name)) {
      ft2 <- get(var_name)
    } else if (is.null(lig_seurat[[input$cellType2]])) {
      ft2 <- fst(paste0("dataFiles/celltype_fst_additionalcompress/200417-ligands-alldata-seurat-p3-",
                        input$cellType2, ".rds.fst", sep = ''))
    }
    plot_df = ft2[c(input$violinInputGene, 'celltype', 'sample', 'rep')]
    plot_df_melt = reshape2::melt(plot_df, id.vars = names(subset(plot_df, select = c('celltype','sample','rep'))))
    data$violin_input_data <- plot_df_melt[c('celltype', 'sample', 'rep', 'variable', 'value')]
    data$violin_input_data <- subset(data$violin_input_data, sample %in% input$violin_inputCytokines)

    violin_plot_results <- lapply(input$violinInputGene, function(myGene) {
      sub <- subset(data$violin_input_data, variable == myGene)
      filtered_df_inp <- sub[sub$sample %in% input$violin_inputCytokines,] # leaving it as is and not deleting
      
      color_palette <- hue_pal()(length(cytokineList))
      
      # Create a ggplot violin plot with scaled width
      p <- ggplot(filtered_df_inp, aes(x = sample, y = value, fill = sample)) +
        geom_violin(scale = "width", size = 0) +
        theme_minimal() +
        theme(axis.text.x = element_text(size = 4.5, angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 6),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size = 7, color = "black"),
              strip.text = element_text(face = "italic", color = "black", size = 5.5),
              strip.placement = "inside",
              strip.background = element_blank(),
              panel.spacing = unit(2, "mm"),
              legend.position = "none") + 
        scale_fill_manual(values = color_palette) +
        facet_grid(variable~., scale = "free_y") +
        scale_y_continuous(breaks = c(0, 3, 6)) +
        ylab("Expression")
      
      fig <- ggplotly(p)

      fig2 <- plot_ly(type = "scatter",
                      mode = "markers",height = 200*length(myGene))%>%
        layout(xaxis=list(visible="FALSE",
                          color="rgba(0,0,0,0)",
                          tickfont =list(color="rgba(0,0,0,0)"),
                          showgrid = F),
               yaxis=list(visible="FALSE",
                          color="rgba(0,0,0,0)",
                          tickfont =list(color="rgba(0,0,0,0)"),
                          showgrid = F),
               showlegend = FALSE)
      fig3 <- subplot(fig, fig2, nrows = 2, heights = c(0.8, 0.2))
      fig3
    })
    violin_data$selectedGenes = input$violinInputGene
    data$violin_plot_input = subplot(violin_plot_results, nrows = length(input$violinInputGene))
    
    showNotification(HTML(paste0(
      'To zoom, select the section you would like to view in detail.')), 
      type = 'message')
    
    # final_data <- list()
    # 
    # for (i in 1:length(input$violinInputGene)) {
    #   
    #   sub <- subset(data$violin_input_data, variable == input$violinInputGene[[i]])
    #   filtered_df_inp <- sub[sub$sample %in% input$violin_inputCytokines,]
    #   color_palette <- hue_pal()(length(cytokineList))
    #   p <- ggplot(filtered_df_inp, aes(x = sample, y = value, fill = sample)) +
    #     geom_violin(scale = "width", size = 0) +
    #     theme_minimal() +
    #     theme(axis.text.x = element_text(size = 4.5, angle = 90, vjust = 0.5, hjust=1),
    #           axis.text.y = element_text(size = 6),
    #           axis.title.x = element_blank(),
    #           axis.title.y = element_text(size = 7, color = "black"),
    #           strip.text = element_text(face = "italic", color = "black", size = 5.5),
    #           strip.placement = "inside",
    #           strip.background = element_blank(),
    #           panel.spacing = unit(2, "mm"),
    #           legend.position = "none") +
    #     scale_fill_manual(values = color_palette) +
    #     facet_grid(variable~., scale = "free_y") +
    #     scale_y_continuous(breaks = c(0, 3, 6)) +
    #     ylab("Expression")
    #   fig <- ggplotly(p)
    #   
    #   fig2 <- plot_ly(height = 200*length(input$violinInputGene))%>%
    #     layout(xaxis=list(visible="FALSE",
    #                       color="rgba(0,0,0,0)",
    #                       tickfont =list(color="rgba(0,0,0,0)"),
    #                       showgrid = F),
    #            yaxis=list(visible="FALSE",
    #                       color="rgba(0,0,0,0)",
    #                       tickfont =list(color="rgba(0,0,0,0)"),
    #                       showgrid = F),
    #            showlegend = FALSE)
    #   
    #   fig3 <- subplot(fig, fig2, nrows = 2, heights = c(0.8, 0.2))
    #   
    #   final_data[[i]] <- fig3
    #   
    #   violin_data$selectedGenes = input$violinInputGene
    # }
    # 
    # data$violin_plot_input = subplot(final_data,
    #                                  nrows = length(input$violinInputGene))
    # 
    # showNotification(HTML(paste0(
    #   'To zoom, select the section you would like to view in detail.')), 
    #   type = 'message')
  }
  else if (input$violin_type_cytokines == 'inp' & length(input$violin_inputCytokines) <= 0){
    showNotification('Please enter in cytokines.', type = 'error')
  }
  else if (input$violin_type_cytokines == 'all'){
    var_name <- paste0("ft2_", input$cellType2)
    if (exists(var_name)) {
      ft2 <- get(var_name)
    } else if (is.null(lig_seurat[[input$cellType2]])) {
      #ft2 <- fst(paste0("dataFiles/celltype_fst_additionalcompress/200417-ligands-alldata-seurat-p3-", 
      ft2 <- fst(paste0("dataFiles/celltype_fst/200417-ligands-alldata-seurat-p3-", 
                                                                          input$cellType2, ".rds.fst", sep = ''))
    }
    plot_df = ft2[c(input$violinInputGene, 'celltype', 'sample', 'rep')]
    plot_df_melt = reshape2::melt(plot_df, id.vars = names(subset(plot_df, select = c('celltype','sample','rep'))))
    data$violin_all_data <- plot_df_melt[c('celltype', 'sample', 'rep', 'variable', 'value')]

    preloaded_figures <- lapply(input$violinInputGene, function(myGene) {
      #fig <- qread(file.path("dataFiles/violin_plot_creation", input$cellType2, paste0(input$cellType2, "_violin_plot_", myGene, ".qs")))
      fig <- qread(file.path("/datadir/violin_plot_creation",input$cellType2, paste0(input$cellType2, "_violin_plot_", myGene, ".qs")))
    }) 
    #preloaded_figures %>% lapply(function(fig) {
      #plotly_data(fig)
    #}) %>% bind_rows() -> data$violin_all_data
    violin_plot_results <- lapply(preloaded_figures, function(fig) {
      # filtered_df_all <- subset(plot_df_melt, variable == myGene)
      # filtered_df_all[filtered_df_all<0] = 0
      # color_palette <- hue_pal()(length(cytokineList))
      # p <- ggplot(filtered_df_all, aes(x = sample, y = value, fill = sample)) +
      #   geom_violin(scale = "width", size = 0) +
      #   theme_minimal() +
      #   theme(axis.text.x = element_text(size = 4.5, angle = 90, vjust = 0.5, hjust=1),
      #         axis.text.y = element_text(size = 6),
      #         axis.title.x = element_blank(),
      #         axis.title.y = element_text(size = 7, color = "black"),
      #         strip.text = element_text(face = "italic", color = "black", size = 5.5),
      #         strip.placement = "inside",
      #         strip.background = element_blank(),
      #         panel.spacing = unit(2, "mm"),
      #         legend.position = "none") +
      #   scale_fill_manual(values = color_palette) +
      #   facet_grid(~variable, scale = "free_y") +
      #   scale_y_continuous(breaks = c(0, 3, 6)) +
      #   ylab("Expression")
      # fig <- ggplotly(p)
      # 
      fig2 <- plot_ly(type = "scatter",
    mode = "markers",height = 200*length(input$violinInputGene))%>%
        layout(xaxis=list(visible="FALSE",
                          color="rgba(0,0,0,0)",
                          tickfont =list(color="rgba(0,0,0,0)"),
                          showgrid = F),
               yaxis=list(visible="FALSE",
                          color="rgba(0,0,0,0)",
                          tickfont =list(color="rgba(0,0,0,0)"),
                          showgrid = F),
               showlegend = FALSE)
      fig3 <- subplot(fig, fig2, nrows = 2, heights = c(0.8, 0.2))
      fig3
    })

    violin_data$selectedGenes = input$violinInputGene
    data$violin_plot_all = subplot(violin_plot_results, nrows = length(input$violinInputGene))

    showNotification(HTML(paste0(
      'To zoom, select the section you would like to view in detail.')), 
      type = 'message')
    
  }
  else{
    return(NULL)
  }
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  print(elapsed_time)}

)
}) # suppress warnings
