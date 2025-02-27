#' Converts the results from CalculateInteractions() into long table format
#'
#' @param results The list returned by CalculateInteractions()
#'
#' @return list of results in long table format
#'
#' @noRd
ConvertToLongTable <- function(results) {
  long_table_list <- NULL
  long_table_names <- names(results)

  for(result in results) {
    colnames(result)[c(4,7,8)] <- c("Gene_secreted","Gene_receptor","pathway_receptor")

    # Get p values and scores for long-table format
    pvalues <- colnames(result)[grepl("_pvalues", colnames(result))]
    scores <- colnames(result)[grepl("_score", colnames(result))]

    # Format the pvalues and scores for each interaction
    interaction_pvalues <- result[ , c("Gene_secreted", "Gene_receptor", "pathway_receptor", pvalues)]
    interaction_scores <- result[ , c("Gene_secreted", "Gene_receptor", "pathway_receptor", scores)]

    # pivot data from wide to long
    df_long <- result %>%
      tidyr::pivot_longer(-c(1:22), names_to = "cell_type_pair", values_to = "value")
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
    colnames(p_value_df)[ncol(p_value_df)] <- "p_val"
    score_df$pval <- p_value_df$p_val
    score_df$receptor <- gsub("_score", "", score_df$receptor)

    long_table_list[[length(long_table_list)+1]] <- score_df
  }

  long_table_list <- setNames(long_table_list, long_table_names)

  return(long_table_list)
}

#' Analysis for multiple samples
#'
#' Analyzes calculated interaction scores for multiple samples by filtering by both percentage expression and the DEG file
#'
#' @param results_list Results from CalculateInteractions()
#' @param DEG_fn Differentially expressed genes file path
#' @param pct_filter Percentage threshold to filter generated Percentage_Expression.csv file for relevant genes
#' @param control_name The name of the control sample, if applicable.
#' @param mutant_name The name of the mutant sample, if applicable.
#'
#' @return Excel File
#'
#' @keywords internal
AnalyzeMultiple <- function(results_list, DEF_fn, pct_filter, control_name, mutant_name) {
  # TODO: This allows for any amount, tho it seems we only need two (control/mutant).
  long_results <- ConvertToLongTable(results_list)

  control_sample <- control_name
  mutant_sample <- mutant_name

  control_interactions <- long_results[[control_sample]] %>% rename(score_control = score, pval_control = pval)
  mutant_interactions <- long_results[[mutant_sample]] %>% rename(score_mutant = score, pval_mutant = pval)

  #FIXME: Look for other cases that might mess up code
  control_interactions$secretor <- gsub("/", "_", control_interactions$secretor)
  control_interactions$receptor <- gsub("/", "_", control_interactions$receptor)
  mutant_interactions$secretor <- gsub("/", "_", mutant_interactions$secretor)
  mutant_interactions$receptor <- gsub("/", "_", mutant_interactions$receptor)

  # DEGS #
  degs <- read.csv(DEF_fn, row.names = 1)

  # Remove LR pairs based off % expression of ligand or receptor
  perc.expr <- read.csv("output/Percentage_Expression.csv")

  ### Filter by % expression ###

  # First, filter perc.expr for relevant genes and pct.1 >= 0.1
  filtered_expr <- perc.expr %>%
    #filter(Homeostasis >= 0.1) %>%
    filter(colnames(perc.expr[2]) >= pct_filter) %>%
    select(Gene, celltype)

  # Separate the filtering for secreted and receptor genes:
  # Get unique combinations of gene and cell_type for secretors and receptors

  secretor_genes <- filtered_expr %>%
    select(Gene, celltype) %>%
    rename(Gene_secreted = Gene, secretor = celltype)

  receptor_genes <- filtered_expr %>%
    select(Gene, celltype) %>%
    rename(Gene_receptor = Gene, receptor = celltype)

  # Perform two left joins to filter interactions:
  filtered_control_interactions <- control_interactions %>%
    inner_join(secretor_genes, by = c("Gene_secreted", "secretor")) %>%
    inner_join(receptor_genes, by = c("Gene_receptor", "receptor"))

  # First, filter perc.expr for relevant genes and pct.2 >= 0.1
  filtered_expr <- perc.expr %>%
    #filter(Recovery_d2 >= 0.1) %>%
    filter(colnames(perc.expr[3]) >= pct_filter) %>%
    select(Gene, celltype)

  # Separate the filtering for secreted and receptor genes:
  # Get unique combinations of gene and cell_type for secretors and receptors

  secretor_genes <- filtered_expr %>%
    select(Gene, celltype) %>%
    rename(Gene_secreted = Gene, secretor = celltype)

  receptor_genes <- filtered_expr %>%
    select(Gene, celltype) %>%
    rename(Gene_receptor = Gene, receptor = celltype)

  # Perform two left joins to filter interactions:
  filtered_mutant_interactions <- mutant_interactions %>%
    inner_join(secretor_genes, by = c("Gene_secreted", "secretor")) %>%
    inner_join(receptor_genes, by = c("Gene_receptor", "receptor"))

  interactions_two_conditions <- filtered_control_interactions %>%
    full_join(filtered_mutant_interactions)

  interactions_two_conditions <- interactions_two_conditions %>%
    mutate(score_control = case_when(is.na(score_control) ~ 0,
                                     .default = score_control)) %>%
    mutate(score_mutant = case_when(is.na(score_mutant) ~ 0,
                                    .default = score_mutant)) %>%
    mutate(pval_control = case_when(is.na(pval_control) ~ 1,
                                    .default = pval_control)) %>%
    mutate(pval_mutant = case_when(is.na(pval_mutant) ~ 1,
                                   .default = pval_mutant))

  ###

  # check which direction DEG
  interactions_two_conditions$secreted_gene_deg_status <- NA
  interactions_two_conditions$receptor_gene_deg_status <- NA

  # loop through secreted genes
  for (i in c(1:nrow(interactions_two_conditions))){
    # Get name of secreted gene and which cell is secreting
    gene_sec <- interactions_two_conditions$Gene_secreted[i]
    sec <- interactions_two_conditions$secretor[i]
    # Assign whether secreted gene from cell is DEG (and which direction is log2fc higher)
    indices <- which(degs$gene == gene_sec & degs$cell_type == sec)

    if (length(indices) == 0){
      next
    }
    if (degs$avg_log2FC[indices] > 0){
      interactions_two_conditions$secreted_gene_deg_status[i] <- mutant_sample
      next
    }
    if (degs$avg_log2FC[indices] < 0){
      interactions_two_conditions$secreted_gene_deg_status[i] <- control_sample
      next
    }
  }

  # loop through receptor genes
  for (i in c(1:nrow(interactions_two_conditions))){
    # Get name of secreted gene and which cell is secreting
    gene_rec <- interactions_two_conditions$Gene_receptor[i]
    rec <- interactions_two_conditions$receptor[i]
    # Assign whether secreted gene from cell is DEG (and which direction is log2fc higher)
    indices <- which(degs$gene == gene_rec & degs$cell_type == rec)

    if (length(indices) == 0){
      next
    }
    if (degs$avg_log2FC[indices] > 0){
      interactions_two_conditions$receptor_gene_deg_status[i] <- mutant_sample
      next
    }
    if (degs$avg_log2FC[indices] < 0){
      interactions_two_conditions$receptor_gene_deg_status[i] <- control_sample
      next
    }
  }

  # Calculate score difference and log2fc
  interactions_two_conditions$score_log2fc <- log2((interactions_two_conditions$score_mutant + 1e-9) / (interactions_two_conditions$score_control + 1e-9)) # use pseudocount of 1e-9
  interactions_two_conditions$score_difference <- interactions_two_conditions$score_mutant - interactions_two_conditions$score_control
  interactions_two_conditions$direction <- ifelse(interactions_two_conditions$score_mutant > interactions_two_conditions$score_control, mutant_sample, control_sample)

  interactions_two_conditions <- interactions_two_conditions %>%
    mutate(specificity = case_when(((pval_mutant < 0.05 & pval_control < 0.05) & ( ! is.na(receptor_gene_deg_status) | ! is.na(secreted_gene_deg_status))) ~ "significant_signaling_both_conditions",
                                   ((pval_mutant < 0.05 & pval_control >= 0.05) & ( ! is.na(receptor_gene_deg_status) | ! is.na(secreted_gene_deg_status))) ~ paste0("significant_signaling_in_", mutant_sample),
                                   ((pval_mutant >= 0.05 & pval_control < 0.05) & ( ! is.na(receptor_gene_deg_status) | ! is.na(secreted_gene_deg_status))) ~ paste0("significant_signaling_in_", control_sample),
                                   ((pval_mutant >= 0.05 & pval_control >= 0.05) & ( ! is.na(receptor_gene_deg_status) | ! is.na(secreted_gene_deg_status))) ~ "tissue_nonspecific_changed_signaling", #TODO What to change name to? Or just remove?
                                   (pval_mutant < 0.05 | pval_control < 0.05) ~ "unchanged_signaling",
                                   .default = "not_significant")
    )

  # Split into 4 data frames based off specificity level
  significant_signaling_both_conditions <- interactions_two_conditions %>% filter(specificity == "significant_signaling_both_conditions")
  significant_signaling_in_Mutant <- interactions_two_conditions %>% filter(specificity == paste0("significant_signaling_in_", mutant_sample))
  significant_signaling_in_Control <- interactions_two_conditions %>% filter(specificity == paste0("significant_signaling_in_", control_sample))
  unchanged_signaling <- interactions_two_conditions %>% filter(specificity == "unchanged_signaling")


  sheet_name <- NULL

  # Create a new workbook
  wb <- createWorkbook()

  # Add sheets with data using set unique names instead of using dynamic expressions
  addWorksheet(wb, sheetName = "Significant_Both_Conditions")
  writeData(wb, sheet = "Significant_Both_Conditions", significant_signaling_both_conditions)

  addWorksheet(wb, sheetName = paste0("Significant_In_", mutant_sample))
  writeData(wb, sheet = paste0("Significant_In_", mutant_sample), significant_signaling_in_Mutant)

  addWorksheet(wb, sheetName = paste0("Significant_In_", control_sample))
  writeData(wb, sheet = paste0("Significant_In_", control_sample), significant_signaling_in_Control)

  addWorksheet(wb, sheetName = "Unchanged")
  writeData(wb, sheet = "Unchanged", unchanged_signaling)

  # Save the workbook
  saveWorkbook(wb, "output/interactions/interactions_multi/results.xlsx", overwrite = TRUE)

  print("AnalyzeMultiple() results saved!")

  list_test <- NULL
  list_test[[length(list_test)+1]] <- significant_signaling_both_conditions
  list_test[[length(list_test)+1]] <- significant_signaling_in_Mutant
  list_test[[length(list_test)+1]] <- significant_signaling_in_Control
  list_test[[length(list_test)+1]] <- unchanged_signaling

  list_test <- setNames(list_test, c("significant_signaling_both_conditions", paste0("significant_signaling_in_", mutant_sample), paste0("significant_signaling_in_", control_sample), "unchanged_signaling"))

  return(list_test)
}
