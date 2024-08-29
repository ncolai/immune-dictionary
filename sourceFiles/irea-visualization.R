
IreaNetworkCircosAll = function(df_irea_network) {

  library(circlize)
  library(colorspace)
  circos.clear()

  # Prepare data for circos plot
  cytokine_spreadsheet = read.xlsx("dataFiles/cytokine_list.xlsx")

  unique_cytokines = unique(df_irea_network$Cytokine)
  unique_cytokines = factor(unique_cytokines, levels = setdiff(cytokine_spreadsheet$Cytokine_DisplayName, "Combo"))
  unique_cytokines = levels(droplevels(unique_cytokines))

    num_cytokines = length(unique_cytokines)
  unique_celltypes_unordered = unique(c(df_irea_network$CelltypeProducer, df_irea_network$CelltypeReceiver))
  num_celltypes = length(unique_celltypes_unordered)

  df_irea_network$CytokineNum = mapvalues(df_irea_network$Cytokine, from = unique_cytokines, to = 1:num_cytokines)

  set.seed(0)
  cytokine_colors = rand_color(length(unique_cytokines))
  names(cytokine_colors) = unique_cytokines

  # Load cell type color list
  celltype_spreadsheet = read.xlsx("dataFiles/celltype_list.xlsx")
  # Arrange cell types to be plotted around the circle based on pre-defined list
  unique_celltypes = unique_celltypes_unordered[na.omit(order(match(unique_celltypes_unordered, celltype_spreadsheet$Celltype_DisplayName)))]

  col_list = lighten(celltype_spreadsheet$Celltype_Color, 0.1)
  names(col_list) = celltype_spreadsheet$Celltype_DisplayName #TODO: change to display name

  # Draw circos plot
  layout(matrix(c(1,1,1,1,1,2), nrow = 1, byrow = TRUE))

  circos.clear()

  par(mar = c(1,1,1,1))
  circos.par(cell.padding = c(0, 0, 0, 0))
  sectors = 1:num_celltypes
  circos.initialize(unique_celltypes, xlim = c(0, num_cytokines))
  circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = col_list[unique_celltypes], track.margin = c(0, 0), bg.border = NA,
               panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2],
                             labels = CELL_META$sector.index,
                             facing = "bending.inside",
                             niceFacing = TRUE,
                             cex = 1.8,
                             adj = c(0.5, -0.3))
               })


  for (ii in 1:nrow(df_irea_network)) {
    df_irea_network_ii = df_irea_network[ii, ]
    cytokine_ii = as.character(df_irea_network_ii$Cytokine)
    cytokine_ii_color = cytokine_colors[cytokine_ii]
    circos.link(df_irea_network_ii$CelltypeProducer, as.numeric(as.vector(df_irea_network_ii$CytokineNum)),
                df_irea_network_ii$CelltypeReceiver, as.numeric(as.vector(df_irea_network_ii$CytokineNum)),
                h = 0.4, col = cytokine_ii_color, directional = 1, w = 0, arr.width = 0.1)
  }

  circos.clear()


  # Plot legend
  par(mar = c(0,1,0,0))
  plot(x=rep(0, num_cytokines), y=1:num_cytokines, col = cytokine_colors[as.character(unique_cytokines)], pch = 15,
       xaxt='n', yaxt='n', ann=FALSE, bty="n", xlim = c(0, 1), ylim = c(num_cytokines, 1))
  text(x=rep(0, num_cytokines), y=1:num_cytokines, labels = as.character(unique_cytokines), col="black", pos = 4)



  par(mfrow = c(1, 1)) # Put plotting arrangement back to original state

}






IreaNetworkCircosIndividual = function(df_irea_network, n_cols = 5, cytokines_to_plot = "all") {

  library(circlize)
  circos.clear()

  # Select the cytokines to plot
  if (cytokines_to_plot[1] != "all") {
    df_irea_network = subset(df_irea_network, Cytokine %in% cytokines_to_plot)
  }

  # Prepare data for circos plot
  cytokine_spreadsheet = read.xlsx("dataFiles/cytokine_list.xlsx")

  unique_cytokines = unique(df_irea_network$Cytokine)
  unique_cytokines = factor(unique_cytokines, levels = setdiff(cytokine_spreadsheet$Cytokine_DisplayName, "Combo"))
  unique_cytokines = levels(droplevels(unique_cytokines))
  num_cytokines = length(unique_cytokines)
  unique_celltypes_unordered = unique(c(df_irea_network$CelltypeProducer, df_irea_network$CelltypeReceiver))
  num_celltypes = length(unique_celltypes_unordered)

  df_irea_network$CytokineNum = mapvalues(df_irea_network$Cytokine, from = unique_cytokines, to = 1:num_cytokines)

  # set.seed(0)
  # cytokine_colors = rand_color(length(unique_cytokines))
  # names(cytokine_colors) = unique_cytokines

  # Load the cell type color list
  celltype_spreadsheet = read.xlsx("dataFiles/celltype_list.xlsx")
  # Arrange cell types to be plotted around the circle based on pre-defined list
  unique_celltypes = unique_celltypes_unordered[na.omit(order(match(unique_celltypes_unordered, celltype_spreadsheet$Celltype_DisplayName)))]

  col_list = lighten(celltype_spreadsheet$Celltype_Color, 0.1)
  names(col_list) = celltype_spreadsheet$Celltype_DisplayName

  # Set the layout for plotting many circos plots
  layout(matrix(1:(ceiling(num_cytokines / n_cols) * n_cols), ncol = n_cols, byrow = TRUE))

  # Make a separate plot for each unique cytokine
  for (cc in unique_cytokines) {
    par(mar = c(0, 0, 1, 0))
    circos.par(cell.padding = c(0, 0, 0, 0))
    sectors = 1:num_celltypes
    circos.initialize(unique_celltypes, xlim = c(0, num_cytokines))
    circos.track(ylim = c(0, 1), track.height = 0.15,
                 bg.col = col_list[unique_celltypes], bg.border = NA)

    df_irea_network_cc = subset(df_irea_network, Cytokine == cc)
    for (ii in 1:nrow(df_irea_network_cc)) {
      df_irea_network_ii = df_irea_network_cc[ii, ]
      circos.link(df_irea_network_ii$CelltypeProducer, as.numeric(as.vector(df_irea_network_ii$CytokineNum)),
                  df_irea_network_ii$CelltypeReceiver, as.numeric(as.vector(df_irea_network_ii$CytokineNum)),
                  h = 0.4, col = "blue", directional = 1, w = 0, arr.width = 0.1, arr.length = 0.2)
    }
    title(cc, font.main = 1)
    circos.clear()
  }

  par(mfrow = c(1, 1)) # Put plotting arrangement back to original state
}



