

# Add to the df_irea whether the cytokine genes are expressed in the specific cell type in the specific sample
AddLigandExpression = function(input_profile, df_irea, threshold = 0.05, species = "mouse") {

  cytokine_receptor_map = read.xlsx("dataFiles/Fig4a-ExtFig22-gene-protein-map.xlsx")
  cytokine_protein_list = unique(df_irea$Cytokine)
  cytokine_protein_list = cytokine_protein_list[!is.na(cytokine_protein_list)]

  if (species == "human") {
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

  # df_irea_tmp = df_irea[, c("Cytokine", "LigandExpressed")]
  # df_irea_tmp = df_irea_tmp[!duplicated(df_irea_tmp), ]

  return(df_irea)

}





AddReceptorExpression = function(input_profile, df_irea, threshold = 0.1, species = "mouse") {

  cytokine_receptor_map = read.xlsx("dataFiles/Fig4a-ExtFig22-gene-protein-map.xlsx")
  cytokine_protein_list = unique(df_irea$Cytokine)
  cytokine_protein_list = cytokine_protein_list[!is.na(cytokine_protein_list)]

  if (species == "human") {
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

  # df_irea_tmp = df_irea[, c("Cytokine", "ReceptorExpressed")]
  # df_irea_tmp = df_irea_tmp[!duplicated(df_irea_tmp), ]

  return(df_irea)
}





