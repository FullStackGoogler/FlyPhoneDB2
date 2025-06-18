#' Raw Interaction Score Calculation
#'
#' Reads in user inputted counts/metadata file pairing and calculates
#' the raw interaction score for all cell-cell combinations.
#'
#' @param counts_fn The file path of the gene/cell matrix containing all expression counts
#' @param metadata_fn The file path of the accompanying metadata file for counts
#' @param LR_pairs The ligand/receptor annotation object chosen based on the knowledgebase version
#' @param pathway_components The pathway core components object chosen based on the knowledgebase version
#' @param perm_times The number of times to shuffle the clusters and calculate ligand/receptor averages for pvalue calculations
#' @param knowledgebase_version The chosen knowledgebase version
#' @param delimitor The separator for the counts file if it is a .txt file. Default is tab.
#' @param seuratObject Counts/Metadata stored in a Seurat Object.
#' @param base_output_dir The directory for FlyPhone to send all result files to. Defaults to the current working directory.
#'
#' @return NULL
#'
#' @keywords internal
CalculateInteractions <- function(counts_fn, metadata_fn, LR_pairs, pathway_components, perm_times, knowledgebase_version, delimitor, seuratObject, base_output_dir) {
  # Loading Data ---------------------------------------------------------------

  output_dir <- paste0(base_output_dir, "output")
  output_base_filename <- "interactions-table-format-unfiltered_"

  # Whether or not we analyse more than one dataset
  multiple_samples <- FALSE

  counts <- NULL
  metadata <- NULL

  counts_split <- NULL
  metadata_split <- NULL

  results <- NULL
  results_names <- NULL

  if(!is.null(seuratObject) && file.exists(seuratObject)) {
    seuratObj <- readRDS(seuratObject)

    metadata <- seuratObj[[]]
    metadata <- rownames_to_column(metadata, "X")
    counts <- seuratObj[["RNA"]]$counts

    if(ncol(metadata) == 2) {
      colnames(metadata)[2] <- "celltype" # This only exists so that the old FlyPhone V1 website example file works
    } else if(ncol(metadata) > 2) {
      metadata <- metadata %>%
        rename("celltype" = cluster)
    }

    if( ("Condition" %in% colnames(metadata)) & (length(unique(metadata$Condition)) > 1) ) {
      print("Splitting metadata...")
      multiple_samples <- TRUE
      metadata_split <- split(metadata, metadata$Condition)
    }
  } else {
    print("Reading in metadata...")

    metadata <- read.csv(metadata_fn)

    if(ncol(metadata) == 2) {
      colnames(metadata)[2] <- "celltype"
    } else if(ncol(metadata) > 2) {
      metadata <- metadata %>%
        rename("celltype" = cluster)
    }

    if( ("Condition" %in% colnames(metadata)) & (length(unique(metadata$Condition)) > 1) ) {
      print("Splitting metadata...")
      multiple_samples <- TRUE
      metadata_split <- split(metadata, metadata$Condition)
    }

    print("Reading in counts...")

    # Read in Counts
    file_type = tools::file_ext(counts_fn)

    if(identical(file_type, "txt")) { # Textfile
      counts <- read.table(counts_fn, header = TRUE, sep = delimitor, row.names = 1, check.names = FALSE)
    } else if(identical(file_type, "csv")) { # CSV file
      counts <- read.csv(counts_fn, row.names = 1, check.names = FALSE)
    } else {
      stop("Unsupported filetype for the count matrix. Please upload either a .CSV or properly delimited .TXT file.")
    }
  }

  print("Creating SeuratObjects...")

  # SeuratObject(s) Creation ---------------------------------------------------
  seurat_object_list <- list()

  if(multiple_samples) {
    print("Splitting up counts matrix...")

    # Split counts up if there are multiple samples
    counts_split <- splitCounts(metadata_split, counts)

    for(sample in names(metadata_split)) {
      seuratObj <- CreateSeuratObject(counts = counts_split[[sample]], meta.data = metadata_split[[sample]], project = sample)
      seurat_object_list <- append(seurat_object_list, seuratObj)
    }
  } else {
    seuratObj <- CreateSeuratObject(counts = counts, meta.data = metadata)
    seurat_object_list <- append(seurat_object_list, seuratObj)
  }

  print("Creating percentage expression file...")

  # Create percentage expression file for later analysis
  percentage_expression <- list()

  for(i in seq_along(seurat_object_list)) {
    seuratObj <- seurat_object_list[[i]]

    count_pe <- as.data.frame(GetAssayData(seuratObj, layer = "counts"))
    count_pe <- as.data.frame(t(count_pe))

    metdata_pe <- seuratObj@meta.data

    combined_pe <- cbind(metdata_pe[, c("celltype")], count_pe)
    colnames(combined_pe)[1] <- "celltype"

    # Get all column names - celltype to get all gene columns
    gene_cols <- setdiff(names(combined_pe), c("celltype"))

    # Group by celltype and calculate percentage expression
    perc_df <- combined_pe %>%
      group_by(celltype) %>%
      summarize(across(all_of(gene_cols), ~ sum(. > 0) / n(), .names = "{.col}"))

    # Convert to long format
    result_df <- perc_df %>%
      pivot_longer(
        cols = all_of(gene_cols),
        names_to = "Gene",
        values_to = gsub("[_/, ]", "-", seuratObj@project.name)
      )

    result_df <- result_df[, c(2, 3, 1)]

    percentage_expression[[i]] <- result_df
  }

  percentage_expression_reduced <- Reduce(function(x, y) merge(x, y, by = c("Gene", "celltype")), percentage_expression)
  pctexpr_final <- percentage_expression_reduced %>%
    select(Gene, everything()) %>%
    relocate(celltype, .after = last_col())

  write.csv(pctexpr_final, paste0(base_output_dir, ".temp/Percentage_Expression.csv"), row.names = FALSE)

  print("Percentage Expression saved to \".temp/Percentage_Expression.csv\"!")

  # Interaction Score Calculation ----------------------------------------------
  for(seuratObj in seurat_object_list) {
    results_names[[length(results_names)+1]] <- seuratObj@project.name

    print(paste0("Calculating interaction scores for: ", seuratObj@project.name))

    # Seurat Object preprocessing & normalization
    seuratObj <- NormalizeData(seuratObj)
    seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
    all_genes <- rownames(seuratObj)
    seuratObj <- ScaleData(seuratObj, features = all_genes)
    seuratObj$cluster <- as.factor(seuratObj$celltype)
    Idents(seuratObj) <- "cluster"

    # Get expression matrix and meta data from seurat object; perform normalization
    exprMat <- as.matrix(GetAssayData(seuratObj, assay = "RNA", layer = "counts"))
    exprMat <- sweep(exprMat, 2, Matrix::colSums(exprMat), FUN = "/") * 10000
    cellInfo <- seuratObj@meta.data

    # Calculate reporter averages ('TF and reporter' included)
    reporter_genes <- pathway_components[pathway_components$role == "reporter", ]
    tf_reporter_genes <- pathway_components[pathway_components$role == "TF and reporter", ]

    common_reporters <- intersect(reporter_genes$gene, row.names(exprMat))
    reporter_genes <- subset(reporter_genes, gene %in% common_reporters)

    common_tf_reporters <- intersect(tf_reporter_genes$gene, row.names(exprMat))
    tf_reporter_genes <- subset(tf_reporter_genes, gene %in% common_tf_reporters)

    # Read in ligand receptor pair data
    names(LR_pairs)[names(LR_pairs) == 'Gene_secreted'] <- "Gene_ligand"
    names(LR_pairs)[names(LR_pairs) == 'Gene_receptor'] <- "Gene_receptor"
    names(LR_pairs)[names(LR_pairs) == 'Pathway_receptor'] <- "Associated_Receptors"

    # Get list of genes in ligand-receptor database and compare to genes detected in samples
    gene_list <- unique(c(LR_pairs$Gene_ligand, LR_pairs$Gene_receptor))
    common_genes <- intersect(gene_list, row.names(exprMat))

    # Subset ligand-receptor pair information to only include connections with genes present in samples
    LR_pairs <- subset(LR_pairs, Gene_ligand %in% common_genes & Gene_receptor %in% common_genes)

    # Get cells as rows and genes as column
    exprMat <- as.matrix(exprMat)
    exprMat <- t(exprMat)

    # Calculate the mean average expression for each gene type (ligand, receptor, reporter, TF and reporter)

    # Subsetting expression matrix cells to keep only possible ligands from LR_pairs database present in expression matrix
    df_Ligand <- exprMat[ , unique(LR_pairs$Gene_ligand)]
    cluster_df_Ligand <- cbind(cellInfo[ , c("cluster"), drop = FALSE], df_Ligand)
    # Average Ligand counts by each cluster
    df_group_by_cluster_Ligand <- cluster_df_Ligand %>%
      group_by(cluster) %>%
      summarise_all(mean) %>%
      as.data.frame()
    # rename rows to cell type/cluster identity and remove column
    row.names(df_group_by_cluster_Ligand) <- df_group_by_cluster_Ligand$cluster
    df_group_by_cluster_Ligand$cluster <- NULL
    # cell type as columns, genes as rows (values are mean ligand count)
    df_group_by_cluster_Ligand <- t(df_group_by_cluster_Ligand)
    # subsetting all cell types(columns) by genes(rows) which are ligands (values are mean expression)
    ligand_avg <- df_group_by_cluster_Ligand[LR_pairs$Gene_ligand, ] %>% as.data.frame()

    # Repeated for receptor genes
    df_Receptor <- exprMat[ , unique(LR_pairs$Gene_receptor)]
    cluster_df_Receptor <- cbind(cellInfo[ , c("cluster"), drop = FALSE], df_Receptor)
    df_group_by_cluster_Receptor <- cluster_df_Receptor %>%
      group_by(cluster) %>%
      summarise_all(mean) %>%
      as.data.frame()
    row.names(df_group_by_cluster_Receptor) <- df_group_by_cluster_Receptor$cluster
    df_group_by_cluster_Receptor$cluster <- NULL
    df_group_by_cluster_Receptor <- t(df_group_by_cluster_Receptor)
    receptor_avg <- df_group_by_cluster_Receptor[LR_pairs$Gene_receptor, ] %>% as.data.frame()

    # Repeated for reporter genes
    df_Reporter <- exprMat[ , colnames(exprMat) %in% common_reporters]
    cluster_df_Reporter <- cbind(cellInfo[ , c("cluster"), drop = FALSE], df_Reporter)
    df_group_by_cluster_reporter <- cluster_df_Reporter %>%
      group_by(cluster) %>%
      summarise_all(mean) %>%
      as.data.frame
    row.names(df_group_by_cluster_reporter) <- df_group_by_cluster_reporter$cluster
    df_group_by_cluster_reporter$cluster <- NULL
    df_group_by_cluster_reporter <- t(df_group_by_cluster_reporter)
    reporter_avg <- df_group_by_cluster_reporter[reporter_genes$gene, ] %>% as.data.frame()

    # Repeated for TF and reporter genes
    df_TF_Reporter <- exprMat[ , colnames(exprMat) %in% common_tf_reporters]
    cluster_df_TF_Reporter <- cbind(cellInfo[ , c("cluster"), drop = FALSE], df_TF_Reporter)
    df_group_by_cluster_TF_Reporter <- cluster_df_TF_Reporter %>%
      group_by(cluster) %>%
      summarise_all(mean) %>%
      as.data.frame
    row.names(df_group_by_cluster_TF_Reporter) <- df_group_by_cluster_TF_Reporter$cluster
    df_group_by_cluster_TF_Reporter$cluster <- NULL
    df_group_by_cluster_TF_Reporter <- t(df_group_by_cluster_TF_Reporter)
    tf_reporter_avg <- df_group_by_cluster_TF_Reporter[tf_reporter_genes$gene, ] %>% as.data.frame()

    # Make a deep copy of the LR pairs
    LR_pairs_one <- LR_pairs # combine

    # Run the permutation test X times sequentially
    num_permutations <- perm_times #TODO: perm_times = 1 breaks calculation of results later
    permutation_results <- map(1:num_permutations, ~ calculate_permuted_avg(cellInfo, exprMat, LR_pairs))

    # Extract ligand and receptor results into separate dataframes
    ligand_results <- map(permutation_results, "ligand_avg")
    receptor_results <- map(permutation_results, "receptor_avg")

    LR_pairs_one <- LR_pairs # combine

    i <- sort(unique(cellInfo$cluster))[1]
    j <- sort(unique(cellInfo$cluster))[2]

    # Calculate cell-cell interaction scores
    for (i in sort(unique(cellInfo$cluster)) ) {

      LR_pairs_combine <- LR_pairs # combine

      for (j in sort(unique(cellInfo$cluster)) ) {
        print(paste0(i, ">", j))
        # Create temporary copy of LR_pairs
        LR_pairs_tmp <- LR_pairs
        # Multiply each cluster's average ligand expression in i by compared cluster's average receptor expression in j, log1p normalized
        LR_pairs_tmp[[paste0(i, ">", j, "_score")]] <- log1p(ligand_avg[[i]]) * log1p(receptor_avg[[j]])

        df <- data.frame(matrix(ncol = perm_times, nrow = nrow(LR_pairs)))
        for (k in 1:perm_times){
          lig_avg_permutate <- ligand_results[[k]]
          rec_avg_permutate <- receptor_results[[k]]
          df[k] <- log1p(lig_avg_permutate[[i]]) * log1p(rec_avg_permutate[[j]])
        }

        LR_pairs_tmp <- cbind(LR_pairs_tmp, df)

        # Tallies up total amount of comparisons where permutated score is higher than the observed score
        if(identical(knowledgebase_version, "Version 1")) {
          LR_pairs_tmp$result <- rowSums(sapply(LR_pairs_tmp[, 13:ncol(LR_pairs_tmp)], function(x) x > LR_pairs_tmp[[paste0(i, ">", j, "_score")]]))
        } else {
          LR_pairs_tmp$result <- rowSums(sapply(LR_pairs_tmp[, 24:ncol(LR_pairs_tmp)], function(x) x > LR_pairs_tmp[[paste0(i, ">", j, "_score")]]))
        }

        # Calculate p score by taking the sum of comparisons where __ and dividing by times of permutation
        LR_pairs_tmp[[paste0(i, ">", j, "_pvalues")]] <- LR_pairs_tmp$result / num_permutations

        # Subset LR_Pairs_tmp to include all rows, the first 12 columns and the p score # CHANGED TO 1:4
        if(identical(knowledgebase_version, "Version 1")) {
          LR_pairs_tmp <- LR_pairs_tmp[ , c(1:12, ncol(LR_pairs_tmp))]
        } else {
          LR_pairs_tmp <- LR_pairs_tmp[ , c(1:23, ncol(LR_pairs_tmp))]
        }

        # Anywhere the score is equal to 0, set pvalue = 1
        LR_pairs_tmp[LR_pairs_tmp[[paste0(i, ">", j, "_score")]] == 0, paste0(i, ">", j, "_pvalues")] <- 1

        # Score table with ligand, receptor, pathway name, interaction score and p value, gets reset after every outer loop (i)
        LR_pairs_combine <- cbind(LR_pairs_combine,
                                  LR_pairs_tmp[ , c(paste0(i, ">", j, "_score"), paste0(i, ">", j, "_pvalues"))]
        )
        # Dataframe containing all interaction scores and p values
        LR_pairs_one <- cbind(LR_pairs_one,
                              LR_pairs_tmp[ , c(paste0(i, ">", j, "_score"), paste0(i, ">", j, "_pvalues"))]
        )
      }
    }

    print("Finished calculating raw interaction values...")

    # Format gene type average expression dataframes
    ligand_avg <- cbind(rownames(ligand_avg), data.frame(ligand_avg, row.names=NULL))
    names(ligand_avg)[names(ligand_avg) == 'rownames(ligand_avg)'] <- "gene"
    ligand_avg <- data.frame(append(ligand_avg, c(gene_type="ligand"), after=1))

    receptor_avg <- cbind(rownames(receptor_avg), data.frame(receptor_avg, row.names=NULL))
    names(receptor_avg)[names(receptor_avg) == 'rownames(receptor_avg)'] <- "gene"
    receptor_avg <- data.frame(append(receptor_avg, c(gene_type="receptor"), after=1))

    reporter_avg <- cbind(rownames(reporter_avg), data.frame(reporter_avg, row.names=NULL))
    names(reporter_avg)[names(reporter_avg) == 'rownames(reporter_avg)'] <- "gene"
    if(nrow(reporter_avg) > 0) { # V1 core components don't have any reporter data, needs to be appended differently
      reporter_avg <- data.frame(append(reporter_avg, c(gene_type="reporter"), after=1))
    } else {
      reporter_avg[1, ] <- NA
      reporter_avg <- data.frame(append(reporter_avg, c(gene_type="reporter"), after=1))
      reporter_avg <- reporter_avg[0, ]
    }

    tf_reporter_avg <- cbind(rownames(tf_reporter_avg), data.frame(tf_reporter_avg, row.names=NULL))
    names(tf_reporter_avg)[names(tf_reporter_avg) == 'rownames(tf_reporter_avg)'] <- "gene"
    if(nrow(tf_reporter_avg) > 0) { # V1 core components don't have any TF and reporter data, needs to be appended differently
      tf_reporter_avg <- data.frame(append(tf_reporter_avg, c(gene_type="TF and reporter"), after=1))
    } else {
      tf_reporter_avg[1, ] <- NA
      tf_reporter_avg <- data.frame(append(tf_reporter_avg, c(gene_type="TF and reporter"), after=1))
      tf_reporter_avg <- tf_reporter_avg[0, ]
    }

    # Combine all averages together
    all_averages <- rbind(ligand_avg, receptor_avg, reporter_avg, tf_reporter_avg)
    all_averages <- data.frame(append(all_averages, c(sample=seuratObj@project.name), after=0))

    cell_columns <- setdiff(names(all_averages), c("sample", "gene", "gene_type"))

    all_averages_filtered <- all_averages # %>%
    # rowwise() %>%
    # mutate(nonzero_proportion = {
    #   expression_vals <- c_across(all_of(cell_columns))
    #   nonzero_count <- sum(expression_vals != 0)
    #   nonzero_count / length(expression_vals)
    # }) %>%
    # filter(nonzero_proportion > 0) %>% # Exclude those with no cell expression
    # select(-nonzero_proportion)

    all_averages_filtered <- data.frame(append(all_averages_filtered, c(Pathway="other"), after=0))

    # Match by gene=gene & gene_type=role between all_averages_filtered=pathway_components
    all_averages_paths <- all_averages_filtered %>%
      left_join(pathway_components, by=c("gene", "gene_type" = "role")) %>%
      mutate(Pathway = ifelse(is.na(pathway), "unknown", pathway)) %>%
      {
        if(identical(knowledgebase_version, "Version 1")) {
          select(., -c("fbgn", "pathway", "GeneID"))
        } else {
          select(., -c("fbgn", "pathway", "version", "source"))
        }
      }

    sampleName <- gsub("[_/, ]", "-", seuratObj@project.name)

    dir.create(paste0(output_dir, "/", sampleName))
    dir.create(paste0(output_dir, "/", sampleName, "/interaction-scores"))

    # Remove "X" from start of columns with numbers at the start
    colnames(all_averages_paths) <- gsub("^X", "", colnames(all_averages_paths))

    write.csv(all_averages_paths, paste0(paste0(output_dir, "/", sampleName), "/pathway-expression-levels_", sampleName, ".csv"), row.names = FALSE)

    cat("FlyPhone mean expression level file saved to:", paste0(paste0(output_dir, "/", sampleName), "/pathway-expression-levels_", sampleName, ".csv"), "\n")

    # Save output
    dir.create(paste0(output_dir), showWarnings = FALSE)
    output_file_path = paste0(paste0(output_dir, "/", sampleName, "/interaction-scores"), "/", output_base_filename, sampleName, ".csv")
    write.csv(LR_pairs_one, file = output_file_path, row.names = FALSE)

    cat("FlyPhone raw interaction score saved to:", output_file_path, "\n")

    # Save data in a list for future functions
    results[[length(results)+1]] <- LR_pairs_one
  }

  # Save sample names for future functions
  formattedNames <- (sapply(results_names, function(x) gsub("[_/, ]", "-", x)))
  write.table(formattedNames, paste0(base_output_dir, "sample_names.txt"), row.names = FALSE, col.names = FALSE, sep = "\n")
  print("Saved sample_names.txt!")

  # Add names onto results and return
  results <- setNames(results, formattedNames)
  return(results)
}

# Helper functions -------------------------------------------------------------

# Splits the counts matrix by Condition if more than one sample is detected
splitCounts <- function(metadata_split, counts) {
  results <- list()

  for(prefix in names(metadata_split)) {
    curr_metadata <- metadata_split[[prefix]]

    barcodes <- curr_metadata[, 1]

    subset_counts <- counts[, colnames(counts) %in% barcodes, drop = FALSE]

    results[[prefix]] <- subset_counts
  }

  return(results)
}

# Define a function to shuffle cluster assignments and calculate ligand/receptor averages
calculate_permuted_avg <- function(cellInfo, exprMat, LR_pairs) {
  # Shuffle the cluster assignments
  shuffled_cellInfo <- cellInfo
  shuffled_cellInfo$cluster <- sample(shuffled_cellInfo$cluster)

  # Recalculate the ligand and receptor averages with shuffled clusters
  cluster_df_Ligand <- cbind(shuffled_cellInfo[, "cluster", drop = FALSE], exprMat[, unique(LR_pairs$Gene_ligand)])
  cluster_df_Receptor <- cbind(shuffled_cellInfo[, "cluster", drop = FALSE], exprMat[, unique(LR_pairs$Gene_receptor)])

  df_group_by_cluster_Ligand <- cluster_df_Ligand %>%
    group_by(cluster) %>%
    summarise_all(mean) %>%
    as.data.frame()

  row.names(df_group_by_cluster_Ligand) <- df_group_by_cluster_Ligand$cluster
  df_group_by_cluster_Ligand$cluster <- NULL
  df_group_by_cluster_Ligand <- t(df_group_by_cluster_Ligand)

  df_group_by_cluster_Receptor <- cluster_df_Receptor %>%
    group_by(cluster) %>%
    summarise_all(mean) %>%
    as.data.frame()

  row.names(df_group_by_cluster_Receptor) <- df_group_by_cluster_Receptor$cluster
  df_group_by_cluster_Receptor$cluster <- NULL
  df_group_by_cluster_Receptor <- t(df_group_by_cluster_Receptor)

  # Subset based on ligand-receptor pairs
  ligand_avg <- df_group_by_cluster_Ligand[LR_pairs$Gene_ligand, ] %>% as.data.frame()
  receptor_avg <- df_group_by_cluster_Receptor[LR_pairs$Gene_receptor, ] %>% as.data.frame()

  # Return ligand and receptor averages as a list
  return(list(ligand_avg = ligand_avg, receptor_avg = receptor_avg))
}
