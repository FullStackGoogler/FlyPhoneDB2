#' Run FlyPhone pipeline
#'
#' Wrapper function for the FlyPhone pipeline.
#'
#' @param counts_fn The filepath to the counts matrix.
#' @param metadata_fn The filepath to the accompanying metadata file.
#' @param DEG The filepath to the DEG file, if provided.
#' @param pct_filter Percent threshold to filter expression of genes by. Default value of 0.1.
#' @param knowledgebase_version Version of knowledgebase to use. Valid values are: "Version 1", "Version 2 All", "Version 2 High", and "Version 2 High/Moderate".
#' @param perm_times The number of times to shuffle cluster alignments and calculate ligand/receptor averages for p-value calculation in raw interaction scoring. Default value of 1000.
#' @param delimitor The separator for the counts file if it is a .txt file. Default is tab.
#' @param control_name The name of the control sample, if applicable.
#' @param mutant_name The name of the mutant sample, if applicable.
#'
#' @return NULL
#'
#' @examples
#' # Afroditi data (GSE218641)
#' RunFlyPhone("GSE218641_filtered_feature_bc_matrix.txt", "GSE218641_GEO_metadata_Afroditi_Petsakou", knowledgebase_version = "Version 2 All")
#' RunFlyPhone(counts_fn = "GSE218641_filtered_feature_bc_matrix.txt", metadata_fn = "GSE218641_GEO_metadata_Afroditi_Petsakou", DEG = "Recovery_vs_Homeostasis_DEG.csv", pct_filter = 0.05, knowledgebase_version = "Version 2 All", perm_times = 100, delimitor = " ")
#'
#' @export
RunFlyPhone <- function(counts_fn = NULL, metadata_fn = NULL, DEG = NULL, pct_filter = 0.1, knowledgebase_version, perm_times = 1000, delimitor = "\t", seuratObject = NULL, control_name = NULL, mutant_name = NULL) {
  # Preprocessing --------------------------------------------------------------

  # Boolean for checking if a DEG file is provided
  DEG_exists <- !is.null(DEG) && file.exists(DEG)

  # Load in annotation and pathway core component data depending on the knowledgebase version chosen
  # annotationObj <- load_object(getAnnotationFile(knowledgebase_version))
  # pathwayObj <- load_object(getPathwayFile(knowledgebase_version))
  annotationObj <- getAnnotationFile(knowledgebase_version)
  pathwayObj <- getPathwayFile(knowledgebase_version)

  # Initial subfolder creation
  dir.create("output")
  dir.create("output/interactions")
  dir.create("output/interactions/interactions_raw")
  dir.create("output/expression_levels")
  dir.create("output/visualizations")
  dir.create("output/visualizations/sample_visualizations")

  # Pipeline -------------------------------------------------------------------

  # Calculate the raw interaction scores
  results <- CalculateInteractions(counts_fn, metadata_fn, annotationObj, pathwayObj, perm_times, knowledgebase_version, delimitor, seuratObject)

  # Boolean for if multiple samples were detected
  isMultiSample <- length(results) > 1

  # Create the rest of the subfolders
  if(isMultiSample & DEG_exists) {
    dir.create("output/interactions/interactions_multi")
    dir.create("output/visualizations/circleplots")
    dir.create("output/visualizations/chord_diagrams")
    dir.create("output/visualizations/interaction_strength")
  } else {
    dir.create("output/interactions/interactions_long")
    dir.create("output/interactions/interactions_long_filtered")
    dir.create("output/visualizations/heatmaps")
    dir.create("output/visualizations/dotplots")
    dir.create("output/visualizations/circleplots")
    dir.create("output/visualizations/chord_diagrams")
    dir.create("output/visualizations/interaction_strength")
  }

  # Only perform multi-sample analysis if a DEG file is provided
  if(isMultiSample & DEG_exists) {
    list_test <- AnalyzeMultiple(results, DEG, pct_filter, control_name, mutant_name) #FIXME: Might not actually use this variable at all
  } else {
    # Single-sample analysis OR Multi-sample analysis WITHOUT a DEG file provided
    AnalyzeSingle(results, knowledgebase_version)
  }

  doMultivis <- isMultiSample && DEG_exists

  # Generate cell-cell communication visualizations
  GenerateVisualizations(counts_fn, metadata_fn, doMultivis, pathwayObj, delimitor, seuratObject)

  # Generate pathway summary visualizations
  PathwayELVisualizations()

  # Remove the percentage expression file from the final output
  file.remove("./output/Percentage_Expression.csv")

  # Create a textfile specifying the type of analysis done
  if(isMultiSample & DEG_exists) {
    write.table("Multi_DEG", "output_type.txt", row.names = FALSE, col.names = FALSE)
  } else if (isMultiSample & !DEG_exists) {
    write.table("Multi", "output_type.txt", row.names = FALSE, col.names = FALSE)
  } else {
    write.table("Single", "output_type.txt", row.names = FALSE, col.names = FALSE)
  }
}

# Helper Functions -------------------------------------------------------------

# Gets the appropriate ligand/receptor pair dataset based on the knowledgebase version chosen
getAnnotationFile <- function(knowledgebase_version) {
  if(identical(knowledgebase_version, "Version 1")) {
    return(FlyPhone:::V1)
  } else if(identical(knowledgebase_version, "Version 2 High")) {
    return(FlyPhone:::V2H)
  } else if(identical(knowledgebase_version, "Version 2 High/Moderate")) {
    return(FlyPhone:::V2HM)
  } else if(identical(knowledgebase_version, "Version 2 All")) {
    return(FlyPhone:::V2A)
  } else {
    stop("Not a valid knowledgebase version! Options include: \"Version 1\",\"Version 2 All\", \"Version 2 High\", or \"Version 2 High/Moderate\".")
  }
}

# Gets the appropriate pathway core component information file based on the knowledgebase version chosen
getPathwayFile <- function(knowledgebase_version) {
  if(identical(knowledgebase_version, "Version 1")) {
    #return("./annotation/Pathway_core_components_2021vs1_clean.txt")
    return(FlyPhone:::pathway_components_v1)
  } else {
    #return("./annotation/Pathway_core_components_Version2-3_final.txt")
    return(FlyPhone:::pathway_components_v2)
  }
}

# load .rda object into a new variable
load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}
