

#' Run IREA on all samples and cell types in the data
#'
#' \code{IreaAll} Compute enrichment score using Wilcoxon rank-sum test
#' and then analyze for receptor and ligand expression across all input
#' samples and cell types
#'
#' @param    input_profile     Gene expression matrix
#' @param    genediff_cutoff   Only include the genes that are differentially expressed above this threshold
#' between cytokine-treated samples and PBS samples to speed up computation
#' @param    species           Animal species being analyzed (mouse/human)
#' @param    threshold_receptor Cutoff to determine whether or not receptor 
#' should be considered expressed
#' @param    threshold_ligand  Cutoff to determine whether or not receptor 
#' should be considered expressed
#' @param    celltypes         Celltypes to run on
#'
#' @export

IreaAll = function(input_profile, genediff_cutoff = 0.25, species = "mouse",
                   threshold_receptor = 0.05,
                   threshold_ligand = 0.05, celltypes) {
  input_colnames = colnames(input_profile)
  samples = unique(gsub("__.*$", "", input_colnames))
  sample_reference = "Control"
  samples = setdiff(samples, sample_reference)
  if (length(celltypes) == 0) { #by default run all celltypes available
    celltypes = unique(gsub("^.*__", "", input_colnames))
  }

  celltype_spreadsheet = read.xlsx("dataFiles/celltype_list.xlsx")

  df_irea = c()

  for (ss in samples) {
    print (paste0("=================== ", ss, " ===================="))
    df_irea_ss = c()

    for (cc in celltypes) {
      print (paste0("------- ", cc, " --------"))
      input_profile_ss_cc = input_profile[, paste0(ss, "__", cc), drop = FALSE]
      input_profile_ss_cc_reltoctrl = input_profile_ss_cc - input_profile[, paste0(sample_reference, "__", cc), drop = FALSE]
      colnames(input_profile_ss_cc_reltoctrl) = gsub("__.*$", "", colnames(input_profile_ss_cc_reltoctrl))

      df_irea_ss_cc = GetEnrichmentScoreProjection(input_profile_ss_cc_reltoctrl,
                                                   input_celltype = cc, genediff_cutoff = genediff_cutoff, species=species)
      df_irea_ss_cc = AddReceptorExpression(df_irea = df_irea_ss_cc,
                                            input_profile = input_profile_ss_cc, threshold = threshold_receptor, species = species)
      df_irea_ss_cc = AddLigandExpression(df_irea = df_irea_ss_cc,
                                          input_profile = input_profile_ss_cc, threshold = threshold_ligand, species = species)
      df_irea_ss = rbind(df_irea_ss, df_irea_ss_cc)
    }

    df_irea = rbind(df_irea, df_irea_ss)
  }

  # Use cell type display name in the data frame to be returned
  df_irea$Celltype = mapvalues(df_irea$Celltype,
                               from = celltype_spreadsheet$Celltype_OriginalName,
                               to = celltype_spreadsheet$Celltype_DisplayName,
                               warn_missing = FALSE)

  return (df_irea)
}










#' Cell-cell interaction network
#'
#' \code{IreaNetwork} First run IreaAll to get the IREA results on all cell 
#' types in all samples, then construct this cell-cell communication network 
#' based on ligand expression, receptor expression, and significant response
#'
#' @param    df_irea_all       The output of IreaAll
#' @param    require_receptor_threshold   Focus on expressed receptors only
#' @param    sig_threshold     Significance required for analysis
#'
#' @export

IreaNetwork = function(df_irea_all, require_receptor_expression = TRUE, sig_threshold = 0.05) {

  samples = as.character(unique(df_irea_all$Sample))
  celltypes = unique(df_irea_all$Celltype)

  df_irea_network_all = c()

  for (ss in samples) {
    df_irea_ss = subset(df_irea_all, Sample == ss)

    df_irea_ss_production = subset(df_irea_ss, LigandExpressed == "Expressed")
    df_irea_ss_production$CelltypeProducer = df_irea_ss_production$Celltype

    df_irea_ss_response = subset(df_irea_ss, padj < sig_threshold) # Only keep cell types with significant responses to a cytokine
    df_irea_ss_response$CelltypeReceiver = df_irea_ss_response$Celltype
    if (require_receptor_expression) {
      # If user specifies that the receptor needs to be expressed to be considered, subset the response data frame further
      # to only include the ones with receptor expressed and significant responses
      df_irea_ss_response = subset(df_irea_ss_response, ReceptorExpressed == "Expressed")
    }

    # Merge the production cell types with receiving cell types to build a cell-cell interaction network for sample ss
    df_irea_ss_network = merge(df_irea_ss_production[, c("Sample", "CelltypeProducer", "Cytokine"), drop = FALSE],
                               df_irea_ss_response[, c("Sample", "CelltypeReceiver", "ES", "padj", "Cytokine", "ReceptorExpressed"), drop = FALSE])

    df_irea_network_all = rbind(df_irea_network_all, df_irea_ss_network)

  }

  return (df_irea_network_all)

}

