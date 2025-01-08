#### Immune response enrichment analysis functions #####

#' IREA analysis for gene set input - rank sum test method
#'
#' \code{GeneSetEnrichmentScore} Perform statistical test (Wilcoxon rank sum test) for cytokine response enrichment when the user input is a list of genes
#'
#' @param    degs              A list of differentially expresssed genes that user would like to investigate
#' @param    input_celltype    Choose from one of the listed cell types that most resembles the input
#' @param    species           Animal species to analyze (mouse/human)
#'
#' @export

GeneSetEnrichmentScore = function(degs, input_celltype, species = "mouse") {
  `%notin%` = Negate(`%in%`)
  library(openxlsx)
  library(plyr)
  library(dplyr)
  set.seed(0)
  
  celltypes = c("ILC", readLines("sourceFiles/lig_seurat_data.txt"))
  if (input_celltype %notin% celltypes) stop(paste0("cell type must be one of the following: ", paste0(celltypes, collapse = ", ")))
  
  # Load reference data
  lig_seurat = readRDS(paste0("dataFiles/ireaData/", input_celltype, ".RDS"))
  
  # Change to official cytokine name
  cytokine_spreadsheet = read.xlsx("dataFiles/irea_cytokine_list.xlsx")
  lig_seurat$sample = mapvalues(lig_seurat$sample,
                                from = cytokine_spreadsheet$Cytokine_OriginalName,
                                to = cytokine_spreadsheet$Cytokine_DisplayName, warn_missing = FALSE)
  
  lig_seurat$Metadata$sample = mapvalues(lig_seurat$Metadata$sample,
                                    from = cytokine_spreadsheet$Cytokine_OriginalName,
                                    to = cytokine_spreadsheet$Cytokine_DisplayName, 
                                    warn_missing = FALSE)
  
  cytokines = setdiff(sort(unique(lig_seurat$Metadata$sample)), "PBS")
  profiles_cc = as.matrix(lig_seurat$RNA_data)
  
  if (tolower(species) == "human") {
    # Load human references
    homolog_table = readRDS("dataFiles/ref_homolog_table.RData")
    rownames(profiles_cc) = mapvalues(rownames(profiles_cc),
                                      from = homolog_table$GeneSymbol1,
                                      to = homolog_table$GeneSymbol2, warn_missing = FALSE)
    
    # For duplciated entries, take the first entry
    profiles_cc = profiles_cc[!duplicated(rownames(profiles_cc)), ]
    
    # Only keep the genes that have corresponding gene symbols between mouse and human
    genes_tokeep = unique(homolog_table$GeneSymbol2)[unique(homolog_table$GeneSymbol2) %in% rownames(profiles_cc)]
    profiles_cc = profiles_cc[genes_tokeep, ]
  }
  
  # only include the DEGs that overlap between user input and reference data
  degs_to_include = degs[degs %in% rownames(profiles_cc)]
  profiles_cc_genes = profiles_cc[degs_to_include, ]
  
  # Calculate enrichment scores
  scores_tmp = aggregate(t(profiles_cc_genes), by = list(lig_seurat$Metadata$sample), mean)
  rownames(scores_tmp) = scores_tmp$Group.1
  scores_tmp$Group.1 = NULL
  scores = apply(scores_tmp, 1, sum)
  scores = scores - scores["PBS"]
  df_irea = data.frame(scores)
  names(df_irea) = "ES"
  df_irea = df_irea[cytokines,, drop = FALSE]
  
  # Assess the significance of enrichment using the Wilcoxon rank sum test between gene set scores on cytokine treated cells
  # and gene set scores on PBS treated cells
  scores_pvals = sapply(cytokines, function(x){
    test_res = wilcox.test(apply(profiles_cc_genes[, rownames(subset(lig_seurat$Metadata, sample == x))], 2, sum),
                           apply(profiles_cc_genes[, rownames(subset(lig_seurat$Metadata, sample == "PBS"))], 2, sum));
    return(test_res$p.value)
  })
  df_irea$pval = as.vector(scores_pvals)
  
  # Perform multiple hypothesis testing correction on all tests
  df_irea$padj = p.adjust(df_irea$pval, method = "fdr")
  
  # Add pseudocount for log transform
  df_irea$nlog10_padj = pmax(0, -log10(df_irea$padj+10e-50))
  
  # Sort results by p-value
  df_irea$Cytokine = rownames(df_irea)
  df_irea$Cytokine = factor(df_irea$Cytokine, levels = df_irea$Cytokine[order(df_irea$nlog10_padj)])
  df_irea$Celltype = input_celltype
  
  return(df_irea)
  
}

#' IREA Polarization analysis for gene set input - rank sum test method
#'
#' \code{PolarizationGeneSetEnrichmentScore} Perform statistical test (Wilcoxon rank sum test) for polarization enrichment when the user input is a list of genes
#'
#' @param    degs              A list of differentially expresssed genes that user would like to investigate
#' @param    input_celltype    Choose from one of the listed cell types that most resembles the input
#' @param    species           Animal species to analyze (mouse/human)
#'
#' @export

PolarizationGeneSetEnrichmentScore = function(degs, input_celltype, species = "mouse") {
  `%notin%` = Negate(`%in%`)
  library(openxlsx)
  library(plyr)
  library(dplyr)
  set.seed(0)
  
  celltypes = c("ILC", readLines("sourceFiles/lig_seurat_data.txt"))
  if (input_celltype %notin% celltypes) stop(paste0("cell type must be one of the following: ", paste0(celltypes, collapse = ", ")))
  
  
  # Load reference data
  filename_cc = paste0("dataFiles/ireaPolarizationData/", input_celltype, ".RDS")
  
  lig_seurat <- readRDS(filename_cc)
  
  profiles_cc = lig_seurat$RNA_data
  metadata_cc = cbind(lig_seurat$Metadata, lig_seurat$Reductions)
  
  # load reference human data
  if (tolower(species) == "human") {
    homolog_table = readRDS("dataFiles/ref_homolog_table.RData")
    rownames(profiles_cc) = mapvalues(rownames(profiles_cc),
                                      from = homolog_table$GeneSymbol1,
                                      to = homolog_table$GeneSymbol2, warn_missing = FALSE)
    
    # For duplciated entries, take the first entry
    profiles_cc = profiles_cc[!duplicated(rownames(profiles_cc)), ]
    
    # Only keep the genes that have corresponding gene symbols between mouse and human
    genes_tokeep = unique(homolog_table$GeneSymbol2)[unique(homolog_table$GeneSymbol2) %in% rownames(profiles_cc)]
    profiles_cc = profiles_cc[genes_tokeep, ]
  }
  
  
  # only include the DEGs that overlap between user input and reference data
  degs_to_include = degs[degs %in% rownames(profiles_cc)]
  profiles_cc_genes = profiles_cc[degs_to_include, ]
  
  
  # Annotate polarization states
  polarization_spreadsheet = read.xlsx("dataFiles/list_map_subcluster_polarization.xlsx")
  metadata_cc$polarization_name = mapvalues(paste0(input_celltype, "_", metadata_cc$subcluster),
                                            from = polarization_spreadsheet$Celltype_subcluster,
                                            to = polarization_spreadsheet$Polarization, warn_missing = FALSE)
  
  
  profiles_cc_agg = aggregate(t(profiles_cc_genes), by = list(metadata_cc$polarization_name), mean)
  rownames(profiles_cc_agg) = profiles_cc_agg$Group.1
  profiles_cc_agg$Group.1 = NULL
  
  ## Reference samples only contain polarized subclusters (i.e. remove the unpolarized clusters, which are marked by S01, S02, etc.)
  reference_samples = setdiff(sort(unique(metadata_cc$polarization_name)),
                              grep(paste0(input_celltype, "_S[0-9][0-9]"), unique(metadata_cc$polarization_name), value = TRUE))
  profiles_cc_agg = profiles_cc_agg[reference_samples, ] # Exclude the "None" polarization state
  
  profile_cc_pbs = profiles_cc_genes[, rownames(subset(metadata_cc, sample == "PBS"))]
  profile_cc_pbs_agg = apply(profile_cc_pbs, 1, mean)
  
  
  
  # Calculate enrichment scores
  scores = apply(profiles_cc_agg, 1, sum)
  scores_pbs = sum(profile_cc_pbs_agg)
  scores = scores - scores_pbs
  df_irea = data.frame(scores)
  names(df_irea) = "ES"
  df_irea = df_irea[reference_samples,, drop = FALSE]
  
  # Assess the significance of enrichment using the Wilcoxon rank sum test between gene set scores on cells in each polarization state
  # and gene set scores on PBS treated cells
  scores_pvals = sapply(reference_samples, function(x){
    test_res = wilcox.test(apply(profiles_cc_genes[, rownames(subset(metadata_cc, polarization_name == x))], 2, sum),
                           apply(profiles_cc_genes[, rownames(subset(metadata_cc, sample == "PBS"))], 2, sum));
    return(test_res$p.value)
  })
  df_irea$pval = as.vector(scores_pvals)
  
  # Perform multiple hypothesis testing correction on all tests
  df_irea$padj = p.adjust(df_irea$pval, method = "fdr")
  
  # Add pseudocount for log transform
  df_irea$nlog10_padj = pmax(0, -log10(df_irea$padj+10e-50))
  
  # Sort results by p-value
  df_irea$Polarization = rownames(df_irea)
  df_irea$Celltype = input_celltype
  
  return(df_irea)
  
}

#### Immune response enrichment analysis functions #####

#' IREA analysis for gene set input - Fisher's test method
#'
#' \code{GeneSetEnrichmentHyperTest} Compute enrichment score using the pathway overrepresentation test,
#' the most commonly used method for pathway analysis which uses the hypergeometric test
#'
#' @param    degs              A list of differentially expresssed genes that user would like to investigate
#' @param    input_celltype    Choose from one of the listed cell types that most resemble the input
#' @param    species           Animal species to analyze (mouse/human)
#'
#' @export

GeneSetEnrichmentHyperTest = function(degs, input_celltype, species = "mouse") {
  `%notin%` = Negate(`%in%`)
  library(openxlsx)
  library(plyr)
  library(dplyr)
  set.seed(0)
  
  celltypes = c("ILC", readLines("sourceFiles/lig_seurat_data.txt"))
  if (input_celltype %notin% celltypes) stop(paste0("cell type must be one of the following: ", paste0(celltypes, collapse = ", ")))
  
  # Load the pre-computed significantly differentially expressed genes for each cytokine in each cell type
  ref_deg_sig_celltype = read.xlsx("dataFiles/SuppTable3_CytokineSignatures.xlsx", sheet = input_celltype)
  ref_expressed_genes = readRDS("dataFiles/ref_expressed_genes_per_celltype.Rda")
  ref_expressed_genes_celltype = subset(ref_expressed_genes, celltype == input_celltype)
  
  # subset into the celltype of interest
  samples = setdiff(read.xlsx("dataFiles/irea_cytokine_list.xlsx")$Cytokine_OriginalName,
                    c("IL2+IL15", "IL2+IFNg"))
  
  # Perform fisher's exact test
  res_pvals = c()
  res_ess = c()
  
  for (ss in samples) {
    markers_ss = subset(ref_deg_sig_celltype, Cytokine_Str == ss)
    
    num_overlap = length(intersect(degs, markers_ss$Gene))
    num_only_user = length(degs) - num_overlap
    num_only_ref = nrow(markers_ss) - num_overlap
    num_neither = length(setdiff(ref_expressed_genes_celltype$genes_expressed,
                                 c(markers_ss$Gene, degs)))
    
    mat_test = matrix(c(num_overlap, num_only_user, num_only_ref, num_neither), nrow = 2)
    
    test_pval = fisher.test(mat_test)$p.value
    test_es = num_overlap / length(degs)
    
    res_pvals = c(res_pvals, test_pval)
    res_ess = c(res_ess, test_es)
    
  }
  
  # Construct a result matrix
  df_irea = data.frame(Cytokine = samples, ES = res_ess, pval = res_pvals)
  
  # Perform multiple hypothesis testing correction
  df_irea$padj = p.adjust(df_irea$pval, method = "fdr")
  
  # Add pseudocount for log transform
  df_irea$nlog10_padj = pmax(0, -log10(df_irea$padj+10e-50))
  df_irea$Celltype = input_celltype
  
  # Order the results by p-value
  df_irea = df_irea[order(df_irea$pval), ]
  
  # Change to official cytokine name
  cytokine_spreadsheet = read.xlsx("dataFiles/irea_cytokine_list.xlsx")
  df_irea$Cytokine = mapvalues(df_irea$Cytokine,
                                from = cytokine_spreadsheet$Cytokine_OriginalName,
                                to = cytokine_spreadsheet$Cytokine_DisplayName, warn_missing = FALSE)
  
  return(df_irea)
}





#' IREA polarization analysis for gene set input - Fisher's test method
#'
#' \code{PolarizationGeneSetEnrichmentHyperTest} Compute enrichment score using the pathway overrepresentation test,
#' the most commonly used method for pathway analysis which uses the hypergeometric test
#'
#' @param    degs              A list of differentially expresssed genes that user would like to investigate
#' @param    input_celltype    Choose from one of the listed cell types that most resemble the input
#' @param    species           Animal species to analyze (mouse/human)
#'
#' @export

PolarizationGeneSetEnrichmentHyperTest = function(degs, input_celltype, species = "mouse") {
  `%notin%` = Negate(`%in%`)
  library(openxlsx)
  library(plyr)
  library(dplyr)
  set.seed(0)
  
  celltypes = c("ILC", readLines("sourceFiles/lig_seurat_data.txt"))
  if (input_celltype %notin% celltypes) stop(paste0("cell type must be one of the following: ", paste0(celltypes, collapse = ", ")))
  
  celltype_spreadsheet = read.xlsx("dataFiles/celltype_list.xlsx")
  input_celltype_display = mapvalues(input_celltype, from = celltype_spreadsheet$Celltype_OriginalName,
                                     to = celltype_spreadsheet$Celltype_DisplayName, warn_missing = FALSE)
  
  print ("Hyper")
  print (input_celltype_display)
  
  # Load the pre-computed significantly differentially expressed genes for each cytokine in each cell type
  ref_deg_sig_celltype = read.xlsx("dataFiles/SuppTable7_Polarization.xlsx", sheet = input_celltype_display)
  ref_expressed_genes = readRDS("dataFiles/ref_expressed_genes_per_celltype.Rda")
  ref_expressed_genes_celltype = subset(ref_expressed_genes, celltype == input_celltype)
  
  # subset into the celltype of interest
  polarization_states = sort(unique(ref_deg_sig_celltype$Polarization))
  
  # Perform fisher's exact test
  res_pvals = c()
  res_ess = c()
  
  for (ss in polarization_states) {
    markers_ss = subset(ref_deg_sig_celltype, Polarization == ss)
    
    num_overlap = length(intersect(degs, markers_ss$Gene))
    num_only_user = length(degs) - num_overlap
    num_only_ref = nrow(markers_ss) - num_overlap
    num_neither = length(setdiff(ref_expressed_genes_celltype$genes_expressed,
                                 c(markers_ss$Gene, degs)))
    
    mat_test = matrix(c(num_overlap, num_only_user, num_only_ref, num_neither), nrow = 2)
    
    test_pval = fisher.test(mat_test)$p.value
    test_es = num_overlap / length(degs)
    
    res_pvals = c(res_pvals, test_pval)
    res_ess = c(res_ess, test_es)
    
  }
  
  # Construct a result matrix
  df_irea = data.frame(Polarization = polarization_states, ES = res_ess, pval = res_pvals)
  
  # Perform multiple hypothesis testing correction
  df_irea$padj = p.adjust(df_irea$pval, method = "fdr")
  
  # Add pseudocount for log transform
  df_irea$nlog10_padj = pmax(0, -log10(df_irea$padj+10e-50))
  df_irea$Celltype = input_celltype
  
  return(df_irea)
}


#### Immune response enrichment score functions #####

#' IREA analysis for transcriptome matrix input
#'
#' \code{GetEnrichmentScoreProjection} Compute enrichment score using the Wilcoxon rank sum test between
#' cosine similarity scores with cytokine-treated samples and cosine similarity scores with control samples
#'
#' @param    input_profile     Gene expression matrix
#' @param    input_celltype    Choose from one of the listed cell types that most resemble the input
#' @param    genediff_cutoff   Only include the genes that are differentially expressed above this threshold
#' between cytokine-treated samples and PBS samples to speed up computation
#' @param    species           Animal species to analyze (mouse/human)
#' @export
GetEnrichmentScoreProjection = function(input_profile, input_celltype, genediff_cutoff = 0.25, species = "mouse") {
  library(openxlsx)
  library(plyr)
  set.seed(0)
  `%notin%` = Negate(`%in%`)
  celltypes = c("ILC", readLines("sourceFiles/lig_seurat_data.txt"))
  
  # Load reference data
  lig_seurat_sub = readRDS(paste0("dataFiles/ireaData/", input_celltype, ".RDS"))
  
  # Change to official cytokine name
  cytokine_spreadsheet = read.xlsx("dataFiles/irea_cytokine_list.xlsx")
  lig_seurat_sub$sample = mapvalues(lig_seurat_sub$sample,
                                    from = cytokine_spreadsheet$Cytokine_OriginalName,
                                    to = cytokine_spreadsheet$Cytokine_DisplayName, 
                                    warn_missing = FALSE)
  
  lig_seurat_sub$Metadata$sample = mapvalues(lig_seurat_sub$Metadata$sample,
                                    from = cytokine_spreadsheet$Cytokine_OriginalName,
                                    to = cytokine_spreadsheet$Cytokine_DisplayName, 
                                    warn_missing = FALSE)
  print("files mapped")
  print(head(lig_seurat_sub$sample))
  
  profiles_cc = as.matrix(lig_seurat_sub$RNA_data)
  #added section to account for human stuff
  if (tolower(species) == "human") {
    #Load human reference data
    homolog_table = readRDS("dataFiles/ref_homolog_table.RData")
    rownames(profiles_cc) = mapvalues(rownames(profiles_cc),
                                      from = homolog_table$GeneSymbol1,
                                      to = homolog_table$GeneSymbol2, warn_missing = FALSE)
    
    # For duplicated entries, take the first entry
    profiles_cc = profiles_cc[!duplicated(rownames(profiles_cc)), ]
    
    # Only keep the genes that have corresponding gene symbols between mouse and human
    genes_tokeep = unique(homolog_table$GeneSymbol2)[unique(homolog_table$GeneSymbol2) %in% rownames(profiles_cc)]
    profiles_cc = profiles_cc[genes_tokeep, ]
  }
  print(profiles_cc)
  metadata_cc = lig_seurat_sub$Metadata
  print(metadata_cc)
  print("cc loaded")
  
  
  homolog_table = readRDS("dataFiles/ref_homolog_table.RData")
  rownames(profiles_cc) = mapvalues(rownames(profiles_cc),
                                    from = homolog_table$GeneSymbol1,
                                    to = homolog_table$GeneSymbol2, warn_missing = FALSE)
  
  # For duplciated entries, take the first entry
  # print(profiles_cc)
  profiles_cc = profiles_cc[!duplicated(rownames(profiles_cc)), ]
  
  
  row_names <- rownames(profiles_cc)
  
  # Sort the row names
  sorted_row_names <- sort(row_names, decreasing = T)
  
  # Print the sorted row names
  print("Sorted row names:")
  print(sorted_row_names)
  
  rownames(input_profile) = mapvalues(rownames(input_profile),
                                           from = homolog_table$GeneSymbol1,
                                           to = homolog_table$GeneSymbol2, warn_missing = FALSE)
  
  # Choose intersection genes
  genes_common = intersect(rownames(input_profile), rownames(profiles_cc))
  print(genes_common)
  
  # Make the final input and final reference by selecting the common genes
  profiles_cc = profiles_cc[genes_common, ]
  # print(profiles_cc)
  input_profile = input_profile[genes_common, , drop = FALSE] # step confusing
  print(head(input_profile))
  
  # Select the genes that are much different between cytokine treatment and PBS
  profiles_cc_agg = aggregate(t(profiles_cc), by = list(metadata_cc$sample), mean)
  rownames(profiles_cc_agg) = profiles_cc_agg$Group.1
  profiles_cc_agg$Group.1 = NULL
  
  # For a gene to be considered different enough from cytokine treatment and control,
  # at least one treatment condition needs to be different from PBS by a certain threshold
  gene_diff = apply(profiles_cc_agg, 2, function(x){max(x-x['PBS'])})
  
  
  # Only use the genes that are significantly differentially expressed to speed up computation
  genes_large_diff = names(gene_diff)[gene_diff > genediff_cutoff]
  profiles_cc = profiles_cc[genes_large_diff, ]
  input_profile = input_profile[genes_large_diff, , drop = FALSE]
  
  dist_input_mat = cbind(input_profile, profiles_cc)
  
  # Project all input vectors onto the reference panel using the cosine similarity metric
  library(philentropy)
  library(reshape2)
  
  coss <- function(x) {crossprod(x)/(sqrt(tcrossprod(colSums(x^2))))} # implement cosine distance using cross product function
  
  projection_scores_raw = coss(dist_input_mat)
  
  # Select the part of the distance matrix that correspond to the results
  projection_scores = projection_scores_raw[1:ncol(input_profile), (ncol(input_profile)+1):ncol(dist_input_mat), drop = FALSE]
  rownames(projection_scores) = colnames(input_profile)
  colnames(projection_scores) = colnames(profiles_cc)
  
  # Format the results matrix
  input_samples = colnames(input_profile)
  reference_samples = setdiff(sort(unique(metadata_cc$sample)), "PBS")
  print(input_samples)
  print("results matrix")
  print(reference_samples)
  mat_pval = matrix(NA, length(input_samples), length(reference_samples),
                    dimnames = list(input_samples, reference_samples))
  df_irea = reshape2::melt(mat_pval)
  df_irea$ES = NA
  names(df_irea) = c("Sample", "Cytokine", "pval", "ES")
  
  # Pre-allocate columns in df_irea
  if (!"pval" %in% colnames(df_irea)) {
    df_irea$pval <- NA
  }
  if (!"ES" %in% colnames(df_irea)) {
    df_irea$ES <- NA
  }
  print("es")
  
  # Compute enrichment score (mean difference between conditions) and p-value
  for (ii in input_samples) {
    for (x in reference_samples) {
      cols_cytokine = which(metadata_cc$sample == x);
      cols_pbs = which(metadata_cc$sample == "PBS");
      
      if (length(cols_cytokine) > 10 && length(cols_pbs) > 3) {
        test_res <- wilcox.test(projection_scores[ii, cols_cytokine], projection_scores[ii, cols_pbs])
        pval <- test_res$p.value
        meandiff <- mean(projection_scores[ii, cols_cytokine], na.rm = TRUE) - mean(projection_scores[ii, cols_pbs], na.rm = TRUE)
      } else {
        pval <- NA
        meandiff <- NA}
      print(pval)
      df_irea[df_irea$Sample == ii & df_irea$Cytokine == x, "pval"] <- pval
      df_irea[df_irea$Sample == ii & df_irea$Cytokine == x, "ES"] <- meandiff
      
    }
  }
  print("es2")
  # Perform multiple hypothesis testing correction
  df_irea$padj = p.adjust(df_irea$pval, method = "fdr")
  df_irea$Celltype = input_celltype
  print("es3")
  # Fill missing data with "not-enriched"
  df_irea$padj[is.na(df_irea$padj)] = 1
  df_irea$ES[is.na(df_irea$ES)] = 0
  
  print("All done?")
  return(df_irea)
}









#### Immune response enrichment score functions #####

#' IREA analysis for transcriptome matrix input
#'
#' \code{PolarizationProjection} Compute enrichment score using the Wilcoxon rank sum test between
#' cosine similarity scores with polarization states and cosine similarity scores with control samples
#'
#' @param    input_profile     Gene expression matrix
#' @param    input_celltype    Choose from one of the listed cell types that most resemble the input
#' @param    genediff_cutoff   Only include the genes that are differentially expressed above this threshold
#' between cytokine-treated samples and PBS samples to speed up computation
#' @param    species           Animal species to analyze (mouse/human)
#' @export

PolarizationProjection = function(input_profile, input_celltype, genediff_cutoff = 0.25, species = "mouse") {
  set.seed(0)
  library(plyr)
  library(reshape2)
  library(openxlsx)
  
  `%notin%` = Negate(`%in%`)
  celltypes = c("ILC", readLines("sourceFiles/lig_seurat_data.txt"))
  if (input_celltype %notin% celltypes) stop(paste0("cell type must be one of the following: ", paste0(celltypes, collapse = ", ")))
  
  
  # Load reference data
  
  filename_cc = paste0("dataFiles/ireaPolarizationData/", input_celltype, ".RDS")
  
  lig_seurat <- readRDS(filename_cc)
  profiles_cc = lig_seurat$RNA_data
  metadata_cc = cbind(lig_seurat$Metadata, lig_seurat$Reductions)
  
  #load reference human data
  if (tolower(species) == "human") {
    homolog_table = readRDS("dataFiles/ref_homolog_table.RData")
    rownames(profiles_cc) = mapvalues(rownames(profiles_cc),
                                      from = homolog_table$GeneSymbol1,
                                      to = homolog_table$GeneSymbol2, warn_missing = FALSE)
    
    # For duplciated entries, take the first entry
    profiles_cc = profiles_cc[!duplicated(rownames(profiles_cc)), ]
    
    # Only keep the genes that have corresponding gene symbols between mouse and human
    genes_tokeep = unique(homolog_table$GeneSymbol2)[unique(homolog_table$GeneSymbol2) %in% rownames(profiles_cc)]
    profiles_cc = profiles_cc[genes_tokeep, ]
  }
  
  # Annotate polarization states
  polarization_spreadsheet = read.xlsx("dataFiles/list_map_subcluster_polarization.xlsx")
  metadata_cc$polarization_name = mapvalues(paste0(input_celltype, "_", metadata_cc$subcluster),
                                            from = polarization_spreadsheet$Celltype_subcluster,
                                            to = polarization_spreadsheet$Polarization, warn_missing = FALSE)
  
  # Choose intersection genes
  genes_common = intersect(rownames(input_profile), rownames(profiles_cc))
  
  # Make the final input and final reference by selecting the common genes
  profiles_cc = profiles_cc[genes_common, ]
  input_profile = input_profile[genes_common, , drop = FALSE]
  
  # Select the genes that are much different between cytokine treatment and PBS
  profiles_cc_agg = aggregate(t(profiles_cc), by = list(metadata_cc$polarization_name), mean)
  rownames(profiles_cc_agg) = profiles_cc_agg$Group.1
  profiles_cc_agg$Group.1 = NULL
  ## Reference samples only contain polarized subclusters (i.e. remove the unpolarized clusters, which are marked by S01, S02, etc.)
  reference_samples = setdiff(sort(unique(metadata_cc$polarization_name)),
                              grep(paste0(input_celltype, "_S[0-9][0-9]"), unique(metadata_cc$polarization_name), value = TRUE))
  profiles_cc_agg = profiles_cc_agg[reference_samples, ] # Exclude the "None" polarization state
  
  profile_cc_pbs = profiles_cc[, rownames(subset(metadata_cc, sample == "PBS"))]
  profile_cc_pbs_agg = apply(profile_cc_pbs, 1, mean)
  
  # For a gene to be considered different enough from cytokine treatment and control,
  # at least one treatment condition needs to be different from PBS by a certain threshold
  gene_diff = apply(rbind(profiles_cc_agg, PBS = profile_cc_pbs_agg), 2, function(x){max(x-x["PBS"])})
  
  # Only use the genes that are significantly differentially expressed to speed up computation
  genes_large_diff = names(gene_diff)[gene_diff>genediff_cutoff]
  profiles_cc = profiles_cc[genes_large_diff, ]
  input_profile = input_profile[genes_large_diff, , drop = FALSE]
  
  dist_input_mat = cbind(input_profile, profiles_cc)
  
  # Project all input vectors onto the reference panel using the cosine similarity metric
  coss <- function(x) {crossprod(x)/(sqrt(tcrossprod(colSums(x^2))))} # implement cosine distance using cross product function
  
  projection_scores_raw = coss(dist_input_mat)
  
  # Select the part of the distance matrix that correspond to the results
  projection_scores = projection_scores_raw[1:ncol(input_profile), (ncol(input_profile)+1):ncol(dist_input_mat), drop = FALSE]
  rownames(projection_scores) = colnames(input_profile)
  colnames(projection_scores) = colnames(profiles_cc)
  
  # Format the results matrix
  input_samples = colnames(input_profile)
  print(input_samples)
  mat_pval = matrix(NA, length(input_samples), length(reference_samples),
                    dimnames = list(input_samples, reference_samples))
  df_irea = reshape2::melt(mat_pval)
  df_irea$ES = NA
  names(df_irea) = c("Sample", "Polarization", "pval", "ES")
  
  # Compute enrichment score (mean difference between conditions) and p-value
  for (ii in input_samples) {
    for (x in reference_samples) {
      cols_cytokine = which(metadata_cc$polarization_name == x);
      cols_pbs = which(metadata_cc$sample == "PBS");
      
      if (length(cols_cytokine) > 10) {
        test_res = wilcox.test(projection_scores[ii, cols_cytokine], projection_scores[ii, cols_pbs])
        pval = test_res$p.value
        meandiff = mean(projection_scores[ii, cols_cytokine]) - mean(projection_scores[ii, cols_pbs])
      } else {
        pval = NULL;
        meandiff = NULL}
      
      df_irea[df_irea$Sample == ii & df_irea$Polarization == x, "pval"] = pval
      df_irea[df_irea$Sample == ii & df_irea$Polarization == x, "ES"] = meandiff
      
    }
  }
  
  # Perform multiple hypothesis testing correction
  df_irea$padj = p.adjust(df_irea$pval, method = "fdr")
  df_irea$Celltype = input_celltype
  
  # Fill missing data with "not-enriched"
  df_irea$padj[is.na(df_irea$padj)] = 1
  df_irea$ES[is.na(df_irea$ES)] = 0
  
  return(df_irea)
}

#' Create ligand expression column
#'
#' \code{AddLigandExpression} Add to the df_irea whether the cytokine ligand 
#' are expressed in the specific cell type in the specific sample
#'
#' @param    input_profile     Gene expression matrix
#' @param    df_irea           IREA dataset to append to
#' @param    threshold         Cutoff to determine expression
#' @param    species           Animal species to analyze (mouse/human)
#' @export
AddLigandExpression = function(input_profile, df_irea, threshold = 0.05, species = "mouse") {
  
  cytokine_receptor_map = read.xlsx("dataFiles/Fig4a-ExtFig22-gene-protein-map.xlsx")
  cytokine_protein_list = unique(df_irea$Cytokine)
  cytokine_protein_list = cytokine_protein_list[!is.na(cytokine_protein_list)]
  
  if (tolower(species) == "human") {
    ligandgene_colname = "Human.gene.symbol"
  } else {
    ligandgene_colname = "Mouse.gene.symbol"
  }
  
  cytokine_gene_list = unlist(strsplit(cytokine_receptor_map[, ligandgene_colname], "(; |, )"))
  cytokine_gene_list = unique(cytokine_gene_list[cytokine_gene_list %in% rownames(input_profile)])
  
  df_irea$LigandExpressed = "NotExpressed"
  
  for (ii in cytokine_protein_list) {
    ligandgene_str = cytokine_receptor_map[cytokine_receptor_map$Cytokine == ii, ligandgene_colname]
    ligandgene_str = ligandgene_str[!is.na(ligandgene_str)]
    all_ligandgenes = unlist(strsplit(ligandgene_str, ", "))
    
    any_ligandgene = FALSE
    for (ll in all_ligandgenes) {
      if (ll %in% rownames(input_profile)) {
        ligandgene_ll_expressed = any(input_profile[ll, ] > threshold) # TODO: change from Control column to the actual sample column
        any_ligandgene = any(any_ligandgene, ligandgene_ll_expressed)
      }
    }
    df_irea[df_irea$Cytokine == ii, "LigandExpressed"] = ifelse(any_ligandgene, "Expressed", "NotExpressed")
    
  }
  
  return(df_irea)
  
}

#' Create receptor expression column
#'
#' \code{AddReceptorExpression} Add to the df_irea whether the cytokine receptor 
#' are expressed in the specific cell type in the specific sample
#'
#' @param    input_profile     Gene expression matrix
#' @param    df_irea           IREA dataset to append to
#' @param    threshold         Cutoff to determine expression
#' @param    species           Animal species to analyze (mouse/human)
#' @export
AddReceptorExpression = function(input_profile, df_irea, threshold = 0.1, species = "mouse") {
  
  cytokine_receptor_map = read.xlsx("dataFiles/Fig4a-ExtFig22-gene-protein-map.xlsx")
  cytokine_protein_list = unique(df_irea$Cytokine)
  cytokine_protein_list = cytokine_protein_list[!is.na(cytokine_protein_list)]
  
  if (tolower(species) == "human") {
    receptor_colname = "Human.main.receptor.gene(s)[Note.1]"
  } else {
    receptor_colname = "Mouse.main.receptor.gene(s)[Note.1]"
  }
  
  receptor_list = unlist(strsplit(cytokine_receptor_map[, receptor_colname], "(; |, )"))
  receptor_list = unique(receptor_list[receptor_list %in% rownames(input_profile)])
  
  df_irea$ReceptorExpressed = "NotExpressed"
  
  for (ii in cytokine_protein_list) {
    receptor_str = cytokine_receptor_map[cytokine_receptor_map$Cytokine == ii, receptor_colname]
    receptor_str = receptor_str[!is.na(receptor_str)]
    all_receptors = unlist(strsplit(receptor_str, "; "))
    
    any_receptor = FALSE
    for (rr in all_receptors) {
      receptor_rr_genes = unlist(strsplit(rr, ", "))
      if (all (receptor_rr_genes %in% rownames(input_profile))) {
        receptor_rr_expressed = all(input_profile[receptor_rr_genes, ] > threshold) #TODO: define the correct column (we are running one column at a time...)
        any_receptor = any(any_receptor, receptor_rr_expressed)
      }
    }
    df_irea[df_irea$Cytokine == ii, "ReceptorExpressed"] = ifelse(any_receptor, "Expressed", "NotExpressed")
    
  }
  
  return(df_irea)
}

#' Extract top genes
#'
#' \code{GetTopGenes} Get the top genes responsible for the 
#' enrichment for the projection method
#'
#' @param    input_profile     Gene expression matrix
#' @param    input_celltype    Choose from one of the listed cell types that most resembles the input
#' @param    cytokine          Specific cytokine to find top genes for
#' @param    genediff_cutoff   Only include the genes that are differentially expressed above this threshold
#' @param    species           Animal species to analyze (mouse/human)
#' @export
GetTopGenes = function(input_profile, input_celltype, cytokine = "IL12", genediff_cutoff = 0.25, species = "mouse") {
  set.seed(0)
  library(plyr)
  library(reshape2)
  library(openxlsx)
  
  `%notin%` = Negate(`%in%`)
  celltypes = c("ILC", readLines("sourceFiles/lig_seurat_data.txt"))
  if (input_celltype %notin% celltypes) stop(paste0("cell type must be one of the following: ", paste0(celltypes, collapse = ", ")))
  
  
  # Load reference data
  filename_cc = paste0("dataFiles/ireaPolarizationData/", input_celltype, ".RDS")
  
  lig_seurat <- readRDS(filename_cc)
  profiles_cc = lig_seurat$RNA_data
  metadata_cc = cbind(lig_seurat$Metadata, lig_seurat$Reductions)
  
  # load human reference data
  if (tolower(species) == "human") {
    homolog_table = readRDS("dataFiles/ref_homolog_table.RData")
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
  
  # Annotate polarization states
  polarization_spreadsheet = read.xlsx("dataFiles/list_polarization_states.xlsx")
  metadata_cc$polarization_name = mapvalues(paste0(input_celltype, "_", metadata_cc$polarization),
                                            from = polarization_spreadsheet$Subcluster_number,
                                            to = polarization_spreadsheet$Polarization_state, warn_missing = FALSE)
  
  # Choose intersection genes
  genes_common = intersect(rownames(input_profile), rownames(profiles_cc))
  
  # Make the final input and final reference by selecting the common genes
  profiles_cc = profiles_cc[genes_common, ]
  input_profile = input_profile[genes_common, , drop = FALSE]
  
  # Select the genes that are much different between cytokine treatment and PBS
  profiles_cc_agg = aggregate(t(profiles_cc), by = list(metadata_cc$polarization_name), mean)
  rownames(profiles_cc_agg) = profiles_cc_agg$Group.1
  profiles_cc_agg$Group.1 = NULL
  reference_samples = setdiff(sort(unique(metadata_cc$polarization_name)), paste0(input_celltype, "_None"))
  profiles_cc_agg = profiles_cc_agg[reference_samples, ] # Exclude the "None" polarization state
  
  profile_cc_pbs = profiles_cc[, rownames(subset(metadata_cc, sample == "PBS"))]
  profile_cc_pbs_agg = apply(profile_cc_pbs, 1, mean)
  
  # For a gene to be considered different enough from cytokine treatment and control,
  # at least one treatment condition needs to be different from PBS by a certain threshold
  gene_diff = apply(rbind(profiles_cc_agg, PBS = profile_cc_pbs_agg), 2, function(x){max(x-x["PBS"])})
  
  # Only use the genes that are significantly differentially expressed to speed up computation
  genes_large_diff = names(gene_diff)[gene_diff>genediff_cutoff]
  profiles_cc = profiles_cc[genes_large_diff, ]
  input_profile = input_profile[genes_large_diff, , drop = FALSE]
}
