#' Analysis for single samples or multiple samples without a DEG file provided
#'
#' Analyzes calculated interaction scores by filtering by significant p-values
#'
#' Converts the results from CalculateInteractions() into two CSV files:
#' - 1: A long table of all interaction scores
#' - 2: A filtered version with interaction scores that have a p-value less than 0.05
#'
#' @param results_list Results from CalculateInteractions()
#'
#' @return CSV Files
#'
#' @keywords internal
AnalyzeSingle <- function(results_list, pct_filter, knowledgebase_version) {
  for(name in names(results_list)) {
    interactions_list <- results_list[[name]]

    perc.expr <- read.csv(".temp/Percentage_Expression.csv")

    filtered_expr <- perc.expr %>%
      filter(colnames(perc.expr[name]) >= pct_filter) %>%
      select(Gene, celltype)

    if(identical(knowledgebase_version, "Version 1")) {
      colnames(interactions_list)[c(3, 6, 7)] <- c("Gene_secreted","Gene_receptor","pathway_receptor")
    } else {
      colnames(interactions_list)[c(4, 7, 8)] <- c("Gene_secreted","Gene_receptor","pathway_receptor")
    }

    # Label clusters that are interesting
    interest_cells <- c()

    # Get p values and scores for long-table format
    pvalues <- colnames(interactions_list)[grepl("_pvalues", colnames(interactions_list))]
    scores <- colnames(interactions_list)[grepl("_score", colnames(interactions_list))]

    # Format the pvalues and scores for each interaction
    interaction_pvalues <- interactions_list[ , c("Gene_secreted", "Gene_receptor", "pathway_receptor", pvalues)]
    interaction_scores <- interactions_list[ , c("Gene_secreted", "Gene_receptor", "pathway_receptor", scores)]

    columns_to_remove <- if(knowledgebase_version == "Version 1") {
      c("geneid_secreted", "FBgn_secreted", "geneid_receptor", "FBgn_receptor",
        "standard", "interaction", "Source", "Additional_info_GO")
    } else {
      c("Pair_ID", "GeneID_secreted", "FBgn_secreted",
        "GeneID_receptor", "FBgn_receptor", "Rank", "Version",
        "Source", "Interaction", "Ligand_annotation", "Ligand_signalP_prediction",
        "Receptor_annotation", "Receptor_TM_prediction", "More_information",
        "Mammalian_ligand.receptor.pair.", "Human_ligand_receptor_.pair.s.",
        "Mouse_ligand_receptor_.pair.s.", "ligand.paralogs.", "receptor.paralogs.")
      # c("Pair_ID", "GeneID_secreted", "FBgn_secreted",
      #   "GeneID_receptor", "FBgn_receptor", "Rank", "Version",
      #   "Source", "Interaction", "Ligand_annotation", "Ligand_signalP_prediction",
      #   "Receptor_annotation", "Receptor_TM_prediction", "More_information",
      #   "Mammalian_ligand-receptor pair?", "Human_ligand_receptor_ pair(s)",
      #   "Mouse_ligand_receptor_ pair(s)", "ligand paralogs?", "receptor paralogs?")
    }

    # Only keep genes with significant percentage expression
    interactions_list <- interactions_list %>%
      inner_join(filtered_expr, by = c("Gene" = "Gene_secreted", "celltype" = "secretor")) %>%
      inner_join(filtered_expr, by = c("Gene" = "Gene_secreted", "celltype" = "receptor"))

    # pivot data from wide to long
    df_long <- interactions_list %>%
      select(-all_of(columns_to_remove)) %>%
      tidyr::pivot_longer(-c(pathway_receptor, Gene_secreted, Gene_receptor), names_to = "cell_type_pair", values_to = "value")
    df_long$cell_type_pair <- gsub("^X", "", df_long$cell_type_pair) # remove X from beginning of cluster column names

    df_long$cell_type_pair <- gsub(".", ">", df_long$cell_type_pair, fixed = TRUE) # cell type A > cell type B : cell type sending signal > (to) receiver
    df_long <- tidyr::separate(df_long, cell_type_pair, sep= ">", into = c("secretor", "receptor"))

    # remove p-value row and add back later
    p_value_df <- df_long %>%
      filter(grepl("pvalue", df_long$receptor))
    score_df <- df_long %>%
      filter(! grepl("pvalue", df_long$receptor))

    # Format long data table
    colnames(score_df)[ncol(score_df)] <- "score"
    colnames(p_value_df)[ncol(score_df)] <- "p_val"
    score_df$pval <- p_value_df$p_val
    score_df$receptor <- gsub("_score", "", score_df$receptor)

    # Filter out non-specific between cell type interactions
    score_df_filtered <- score_df %>%
      filter(pval < 0.05)

    write.csv(score_df_filtered, paste0("output/interactions/interactions_long_filtered/interaction_long_filtered_", name, ".csv"), row.names = FALSE)
    write.csv(score_df, paste0("output/interactions/interactions_long/interaction_long_", name,".csv"), row.names = FALSE)

    print("AnalyzeSingle() results saved!")
  }
}
