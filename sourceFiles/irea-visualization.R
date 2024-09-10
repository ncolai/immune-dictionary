#### Immune response enrichment analysis visualization functions #####

#' Compass plot
#'
#' \code{IreaCompassPlot} Visualize IREA results on a compass plot
#'
#' @param    df_irea        IREA results data frame returned from one of the "GeneSetEnrichmentScore",
#' "GeneSetEnrichmentHyperTest", or "GetEnrichmentScoreProjection" methods
#' @param    color_by       Color the IREA compass plot by p-value or effect size; default is p-value
#' (TODO: add support for ES)
#'

IreaCompassPlot = function(df_irea, color_by = c("pval", "ES"), plot_receptor = FALSE) {
  `%notin%` = Negate(`%in%`)
  
  library(ggplot2)
  library(plyr)
  library(ggnewscale)
  
  # TODO: double check this if statement works
  if (plot_receptor == TRUE && ! "ReceptorExpressed" %in% colnames(df_irea)) {
    warning("To plot receptor expression, please run the AddReceptorExpression function first.
            Plotting without receptor expression.")
    plot_receptor = FALSE
  }
  
  
  df_irea$NES = as.numeric(as.vector(df_irea$ES)) / abs(max(as.numeric(as.vector(df_irea$ES))))
  df_irea$nlog10p = -log10(df_irea$padj + 0.00001)
  df_irea$nlog10p[df_irea$nlog10p<0]=0.000000001
  
  # For visualization purposes, visualize positive enrichment with red, negative with blue. Darker colors represent stronger p-value
  df_irea$nlog10p_signed = df_irea$nlog10p * sign(df_irea$NES)
  max_barlimit = max(df_irea$NES)
  min_barlimit = min(df_irea$NES)
  text_angle = c(seq(0,-180,length.out = 43), seq(360, 180, length.out = 43))
  
  # Order cytokines to be displayed by enrichment score
  #unique_cytokines <- unique(df_irea$Cytokine)
  #df_irea$Cytokine = factor(unique_cytokines, levels = unique_cytokines[order(df_irea$NES)])
  df_irea$Cytokine = factor(df_irea$Cytokine, levels = unique(df_irea$Cytokine[order(df_irea$NES)]))
  
  
  # Annotation data frame, which is used to add text (cytokine names) with a specific angle
  text_angle = c(seq(0,-180,length.out = 43), seq(360, 180, length.out = 43))
  text_hjust = c(rep(0, 43), rep(1, 43))
  df_annotate = data.frame(xx = levels(df_irea$Cytokine), yy = max_barlimit*1.1, aa = text_angle + 90, hh = text_hjust,
                           ff = "black")
  
  
  # Make compass plot
  if (plot_receptor) {
    p1 = ggplot() +
      theme_void() +
      geom_segment(data = df_irea, aes(x = Cytokine, xend = Cytokine), y = 0, yend = max_barlimit,
                   arrow = arrow(length = unit(0.12,"cm")), color = "gray30") +
      geom_bar(data = df_irea, aes(x = Cytokine, y = NES, fill = nlog10p_signed), stat = "identity", size = 0.2) +
      geom_text(data = df_annotate, y = max_barlimit*1.1, aes(x = xx, label = xx, angle = aa, hjust = hh), size = 3.2) +
      scale_fill_gradient2(name = expression(-log[10]*italic(P)),
                           mid = "white", high = "red", low = "blue", midpoint = 0, limits = c(-5, 5),
                           breaks = c(-5.0, -2.5, 0, 2.5, 5.0), labels = c("5.0+", "2.5", "0", "2.5", "5.0+")) +
      new_scale_fill() +
      geom_tile(data = df_irea, aes(x = Cytokine, y = max_barlimit*1.05, fill = ReceptorExpressed), color = "black", size = 0.5, width = 0.8, height = 0.05) +
      scale_fill_manual(guide = "none", name = "Receptor\nexpression", values = c("Expressed" = "gray40", "NotExpressed" = "white")) +
      theme(axis.text.x = element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.key.width = unit(10, "pt"),
            legend.key.height = unit(14, "pt"),
            legend.position = c(0.94, 0.25)) +
      # ylim(-0.01, 0.15) + #
      ylim(min_barlimit, (max_barlimit*1.4)) +
      coord_polar()
    
  } else {
    
    p1 = ggplot() +
      theme_void() +
      geom_segment(data = df_irea, aes(x = Cytokine, xend = Cytokine), y = 0, yend = max_barlimit,
                   arrow = arrow(length = unit(0.12,"cm")), color = "gray30") +
      geom_bar(data = df_irea, aes(x = Cytokine, y = NES, fill = nlog10p_signed), stat = "identity", size = 0.2) +
      geom_text(data = df_annotate, y = max_barlimit*1.05, aes(x = xx, label = xx, angle = aa, hjust = hh), size = 3.2) +
      scale_fill_gradient2(name = expression(-log[10]*italic(P)),
                           mid = "white", high = "red", low = "blue", midpoint = 0, limits = c(-5, 5),
                           breaks = c(-5.0, -2.5, 0, 2.5, 5.0), labels = c("5.0+", "2.5", "0", "2.5", "5.0+")) +
      theme(axis.text.x = element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.position = c(1.15, 0.2),
            legend.key.width = unit(10, "pt"),
            legend.key.height = unit(14, "pt"),
            plot.margin = margin(50,50,50,50, unit = "pt"),
            legend.margin = margin(0,0,0,0),
            legend.box.margin = margin(0,0,0,0)) +
      # ylim(-0.01, 0.15) + #
      ylim(min_barlimit, (max_barlimit*1.2)) +
      coord_polar()
  }
  return (p1)
  
}








#' Radar plot
#'
#' \code{IreaRadarPlot} Visualize IREA results on a radar plot. Default visualization for immune cell polarization.
#'
#' @param    df_irea        IREA results data frame returned from one of the "GeneSetEnrichmentScore",
#' "GeneSetEnrichmentHyperTest", or "GetEnrichmentScoreProjection" methods
#' @param    color_by       Color the IREA compass plot by p-value or effect size; default is p-value
#' (TODO: add support for ES)
#' @param    celltype       Cell type analyzed
#' @param    fill_color     Color to be filled on the radar plot
#' @param    analysis_type  Visualize enrichment of immune cell polarization states or individual cytokines

IreaRadarPlot = function(df_irea,
                         color_by = "ES",
                         analysis_type = "Polarization",
                         input_celltype = "Macrophage",
                         fill_color = "default") {
  `%notin%` = Negate(`%in%`)
  
  library(ggplot2)
  library(plyr)
  library(openxlsx)
  
  #df_irea$Polarization = axes_name
  df_irea$nlog10p = -log10(df_irea$padj + 0.00001)
  df_irea$nlog10p[df_irea$nlog10p<0]=0.000000001
  
  # Load the corresponding files from the immune dictionary
  irea_polarization_states = read.xlsx("dataFiles/list_map_subcluster_polarization.xlsx")
  irea_polarization_states_cc = subset(irea_polarization_states, Celltype == input_celltype)
  irea_polarization_states_cc = irea_polarization_states_cc[!duplicated(irea_polarization_states_cc), ]
  
  celltype_spreadsheet = read.xlsx("dataFiles/celltype_list.xlsx")
  celltype_display = celltype_spreadsheet[celltype_spreadsheet$Celltype_OriginalName == input_celltype, "Celltype_DisplayName"][1]
  
  # Load the default cell type color
  if (fill_color == "default") {
    fill_color = celltype_spreadsheet[celltype_spreadsheet$Celltype_OriginalName == input_celltype, "Celltype_Color"][1]
  }
  
  # Process the axes to be plotted
  axes_name = unique(irea_polarization_states_cc$Polarization)
  
  #rownames(df_irea) = df_irea[[analysis_type]]
  
  # OR #
  
  df_irea$UniqueRowName <- paste0(df_irea[, analysis_type], seq_along(df_irea[, analysis_type]))
  rownames(df_irea) <- df_irea$UniqueRowName
  
  #axes_enrichment = df_irea[axes_name, "ES"]
  axes_enrichment = df_irea[df_irea$UniqueRowName, "ES"]
  axes_enrichment[is.na(axes_enrichment)] = 0 # The states not found have an enrichment score of 0
  axes_enrichment[axes_enrichment<0] = 0 # Do not plot negative values
  
  # Rescale the axis from 0.2 to 1 for visualization purposes
  axes_enrichment = (axes_enrichment - min(axes_enrichment))/ (max(axes_enrichment) - min(axes_enrichment)) + 0.3 # 0.3 is 0 on the plot
  # TODO: right now the maximum is 1, but the maximum should be whatever is in the data
  
  # If no polarization state is significant, show an enrichment score of 0 for every state
  if (min(df_irea$padj) > 0.05) {
    axes_enrichment = rep(0.3, length(axes_enrichment)) # 0.3 is 0 on the plot
  }
  
  # Visualization code
  axes_enrichment = c(axes_enrichment, axes_enrichment[1])
  n_axes = length(axes_name)
  c_degrees = seq(0, 2*pi, length = c(n_axes+1)) + 0.5*pi
  df_grid = data.frame(dd = c_degrees,
                       xx = c(cos(c_degrees)*1.3, cos(c_degrees)*0.8, cos(c_degrees)*0.3),
                       yy = c(sin(c_degrees)*1.3, sin(c_degrees)*0.8, sin(c_degrees)*0.3),
                       gg = c(rep("outer", n_axes+1), rep("mid", n_axes+1), rep("inner", n_axes+1)))
  df_radial = data.frame(xx = c(cos(c_degrees[1:n_axes])*0.3, cos(c_degrees[1:n_axes])*1.3),
                         yy = c(sin(c_degrees[1:n_axes])*0.3, sin(c_degrees[1:n_axes])*1.3),
                         gg = c(1:n_axes, 1:n_axes))
  
  df_axis_label = data.frame(xx = c(0.1, 0.1, 0.1),
                             yy = c(0.3, 0.8, 1.3),
                             text_label = c("0", "0.5", "1"))
  
  df_text = data.frame(xx = c(0, cos(c_degrees[1:n_axes])*1.6),
                       yy = c(0, sin(c_degrees[1:n_axes])*1.6),
                       state_text = c(celltype_display, axes_name))
  # Move text labels on the left and right size a little bit further from the plot
  #df_text$xx[df_text$yy < 1.6 & df_text$yy > - 1.6] = df_text$xx[df_text$yy < 1.6 & df_text$yy > - 1.6]*1.1
  
  df_polygon = data.frame(xx = axes_enrichment * cos(c_degrees),
                          yy = axes_enrichment * sin(c_degrees))
  
  p1 = ggplot() +
    theme_void() +
    geom_path(data = df_grid, aes(x = xx, y = yy, group = gg), color = "blue", alpha = 0.5) +
    geom_path(data = df_radial, aes(x = xx, y = yy, group = gg), color = "black") +
    geom_text(data = df_axis_label, aes(x = xx, y = yy, label = text_label), color = "gray50", size = 5) +
    geom_polygon(data = df_polygon, aes(x = xx, y = yy), alpha = 0.8,
                 color = fill_color, fill = fill_color, size = 1) +
    geom_point(data = df_polygon, aes(x = xx, y = yy), color = fill_color, size = 3) +
    geom_text(data = df_text, aes(x = xx, y = yy, label = state_text), hjust = 0.5, size = 5) +
    xlim(-2, 2) +
    ylim(-2, 2) +
    theme(aspect.ratio = 1)
  
  return (p1)
  
}






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
  text(x=rep(0, num_cytokines), y=1:num_cytokines, labels = as.character(unique_cytokines), col="black", pos = 4, cex=1.2)
  
  
  
  par(mfrow = c(1, 1)) # Put plotting arrangement back to original state
  
  return(recordPlot())
}






IreaNetworkCircosIndividual = function(df_irea_network, n_cols = 5, cytokines_to_plot = "all") {
  
  library(circlize)
  circos.clear()
  
  # Select the cytokines to plot
  if (cytokines_to_plot[1] != "all") {
    df_irea_network = subset(df_irea_network, Cytokine %in% cytokines_to_plot)
  }
  
  # Prepare data for circos plot
  #cytokine_spreadsheet = read.xlsx("~/Dropbox/Hacohen/ligands/spreadsheets/cytokine_list.xlsx")
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
  
  return(recordPlot())
}






IreaTopGenes = function(input_profile,
                        input_celltype = "NK_cell",
                        input_cytokine = "IL12",
                        genediff_cutoff = 0.25,
                        species = "mouse",
                        ngenes_display = 30,
                        user_data_str = "User data\nsignature",
                        return_data = FALSE){
  library(openxlsx)
  library(plyr)
  library(reshape2)
  library(ggplot2)
  library(gridExtra)
  
  set.seed(0)
  `%notin%` = Negate(`%in%`)
  celltypes = c("B_cell", "cDC1", "cDC2", "eTAC", "ILC", "Macrophage", "MigDC",
                "Monocyte", "Neutrophil", "NK_cell", "pDC", "T_cell_CD4", "T_cell_CD8",
                "T_cell_gd", "Treg")
  if (input_celltype %notin% celltypes) stop(paste0("cell type must be one of the following: ", paste0(celltypes, collapse = ", ")))
  
  
  # Load reference data
  #lig_seurat = readRDS("../data/lig_seurat.Rda")
  lig_seurat = readRDS("~/Dropbox/RPackages/IREA//data/lig_seurat.Rda")
  lig_seurat_sub = subset(lig_seurat, celltype == input_celltype)
  
  # Change to official cytokine name
  cytokine_spreadsheet = read.xlsx("~/Dropbox/Hacohen/ligands/spreadsheets/cytokine_list.xlsx")
  input_cytokine = mapvalues(input_cytokine, from = cytokine_spreadsheet$Cytokine_OriginalName,
                             to = cytokine_spreadsheet$Cytokine_DisplayName, warn_missing = FALSE)
  lig_seurat_sub$sample = mapvalues(lig_seurat_sub$sample,
                                    from = cytokine_spreadsheet$Cytokine_OriginalName,
                                    to = cytokine_spreadsheet$Cytokine_DisplayName, warn_missing = FALSE)
  
  profiles_cc = as.matrix(lig_seurat_sub@assays[['RNA']]@data)
  
  # TODO: is there a better way to translate gene names from mouse to human?
  if (species == "human") {
    # TODO: can actually save these datasets in RDS as human references to speed up
    homolog_table = readRDS("~/Dropbox/RPackages/IREA//data/ref_homolog_table.RData")
    rownames(profiles_cc) = mapvalues(rownames(profiles_cc),
                                      from = homolog_table$GeneSymbol1,
                                      to = homolog_table$GeneSymbol2, warn_missing = FALSE)
    
    # For duplciated entries, take the first entry
    profiles_cc = profiles_cc[!duplicated(rownames(profiles_cc)), ]
    # Alternative: If there are duplicated entries, take the mean of duplicated entries
    #profiles_cc = aggregate(profiles_cc, by = list(rownames(profiles_cc)), mean)
    
    # Only keep the genes that have corresponding gene symbols between mouse and human
    genes_tokeep = unique(homolog_table$GeneSymbol2)[unique(homolog_table$GeneSymbol2) %in% rownames(profiles_cc)]
    profiles_cc = profiles_cc[genes_tokeep, ]
  }
  
  
  # Choose intersection genes
  genes_common = intersect(rownames(input_profile), rownames(profiles_cc))
  
  # Make the final input and final reference by selecting the common genes
  profiles_cc = profiles_cc[genes_common, ]
  input_profile = input_profile[genes_common, , drop = FALSE]
  
  # Select the genes that are much different between cytokine treatment and PBS
  profiles_cc_agg = aggregate(t(profiles_cc), by = list(lig_seurat_sub@meta.data$sample), mean)
  rownames(profiles_cc_agg) = profiles_cc_agg$Group.1
  profiles_cc_agg$Group.1 = NULL
  
  # For a gene to be considered different enough from cytokine treatment and control,
  # at least one treatment condition needs to be different from PBS by a certain threshold
  gene_diff = apply(profiles_cc_agg, 2, function(x){max(x - x['PBS'])})
  
  profiles_cc_agg_diff = t(apply(profiles_cc_agg, 2, function(x){x - x['PBS']}))
  
  # Only use the genes that are significantly differentially expressed to speed up computation
  genes_large_diff = names(gene_diff)[gene_diff > genediff_cutoff]
  profiles_cc_agg_diff_inputcytokine = profiles_cc_agg_diff[genes_large_diff, input_cytokine, drop = FALSE]
  input_profile = input_profile[genes_large_diff, , drop = FALSE]
  
  product_input_reference = profiles_cc_agg_diff_inputcytokine * input_profile
  
  product_input_reference_sorted = product_input_reference[order(product_input_reference, decreasing = TRUE),, drop = FALSE]
  
  
  # Get the top genes
  top_genes = rownames(head(product_input_reference_sorted, ngenes_display))
  
  # Create a data frame with the top genes and magnitude and direction of signature
  df_topgenes = data.frame(Reference = profiles_cc_agg_diff_inputcytokine[top_genes,],
                           Userdata = input_profile[top_genes,]*2)
  df_topgenes$Gene = rownames(df_topgenes)
  
  df_topgenes_melt = melt(df_topgenes, id.vars = "Gene")
  df_topgenes_melt$Gene = factor(df_topgenes_melt$Gene, levels = rev(df_topgenes$Gene))
  df_topgenes_melt$AbsValue = abs(df_topgenes_melt$value)
  df_topgenes_melt$Sign = sign(df_topgenes_melt$value)
  df_topgenes_melt$Regulation = mapvalues(df_topgenes_melt$Sign, from = as.factor(c(-1, 1)), to = c("Downregulated", "Upregulated"))
  df_topgenes_melt$AbsValueLeftRight = df_topgenes_melt$AbsValue
  df_topgenes_melt[df_topgenes_melt$variable=="Reference", "AbsValueLeftRight"] = -df_topgenes_melt[df_topgenes_melt$variable=="Reference", "AbsValueLeftRight"]
  
  
  # Make a diverging bar plot by creating two separate bar plots, one pointing to the left, one pointint to the right, then merge
  df_topgenes_melt_left = subset(df_topgenes_melt, variable == "Reference")
  df_topgenes_melt_right = subset(df_topgenes_melt, variable == "Userdata")
  
  p_left = ggplot(df_topgenes_melt_left, aes(x = Gene, y = AbsValueLeftRight)) +
    geom_bar(stat = "identity", aes(fill = Regulation)) +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0)) +
    scale_fill_manual(values = c("Downregulated" = "blue", "Upregulated" = "red")) +
    scale_y_continuous(expand = c(0,0)) +
    ylab("Dictionary\nsignature")
  
  center_plot = ggplot(df_topgenes_melt_left, aes(x=Gene, y=AbsValueLeftRight)) +
    geom_point(aes(y=0), size=0, color="white") +
    geom_text(aes(label=Gene, y=0), hjust=0.5, vjust=0.4, fontface = "italic", size = 3.2) +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    ylab("\n")
  
  p_right = ggplot(df_topgenes_melt_right, aes(x = Gene, y = AbsValueLeftRight)) +
    geom_bar(stat = "identity", aes(fill = Regulation)) +
    coord_flip() +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0)) +
    scale_fill_manual(values = c("Downregulated" = "blue", "Upregulated" = "red")) +
    scale_y_continuous(expand = c(0,0)) +
    ylab(user_data_str)
  
  p_out = arrangeGrob(p_left, center_plot, p_right, ncol=3, widths=c(4, 2, 4))
  
  if (return_data) {
    return(list(p_out, df_topgenes_melt))
  }
}
