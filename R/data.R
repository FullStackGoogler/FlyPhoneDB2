#' FlyPhone Ligand-Receptor Pairs Version 1 - 2021
#'
#' Dataset of all ligand-receptor pairs for the major fly signaling pathways.
#'
#' @format ## `Version 1`
#' A data frame with 196 rows and 11 columns:
#'
#' @docType data
#' @source Manually curated
"V1"

#' FlyPhone Ligand-Receptor Pairs Version 2.5 - 2025
#'
#' Dataset of all ligand-receptor pairs for the major fly signaling pathways.
#'
#' @format ## `Version 2 All`
#' A data frame with 1,748 rows and 22 columns:
#'
#' @docType data
#' @source Manually curated
"V2A"

#' FlyPhone Ligand-Receptor Pairs Version 2.5 - 2025
#'
#' Dataset of ligand-receptor pairs for the major fly signaling pathways that have a rank of "High".
#'
#' @format ## `Version 2 High`
#' A data frame with 259 rows and 22 columns:
#'
#' @docType data
#' @source Manually curated
"V2H"

#' FlyPhone Ligand-Receptor Pairs Version 2.5 - 2025
#'
#' Dataset of all ligand-receptor pairs for the major fly signaling pathways.
#'
#' @format ## `Version 2 High/Moderate`
#' A data frame with 910 rows and 22 columns:
#'
#' @docType data
#' @source Manually curated
"V2HM"


#' Pathway Core Components Version 1 - 2021
#'
#' Dataset of core components information for each signaling pathway.
#'
#' @format ## `Pathway Components V1`
#' A data frame with 303 rows and 5 columns:
#'
#' @docType data
#' @source Manually curated
"pathway_components_v1"

#' Pathway Core Components Version 2.5 - 2025
#'
#' Dataset of core components information for each signaling pathway.
#'
#' @format ## `Pathway Components V2`
#' A data frame with 458 rows and 6 columns:
#'
#' @docType data
#' @source Manually curated
"pathway_components_v2"


#' Sample data: example_seurat_obj
#'
#' A down-sampled Seurat object (Day 5 Yki snRNA-seq data) for demonstrating package functionality.
#'
#' @format ## `Seurat Object`
#' \describe{
#'   \item{features}{16,443 genes}
#'   \item{samples}{4,187 cells}
#'   \item{assays}{\strong{RNA} (raw counts and normalized data)}
#'   \item{layers}{\strong{counts}: unnormalized expression counts; \strong{data}: log-normalized expression}
#'   \item{reductions}{\strong{umap}: UMAP dimensional reduction calculated}
#' }
#' @source A Seurat object generated from Day 5 Yki snRNA-seq data using 10 randomly selected clusters downsampled to 500 cells in each cluster
#'
#' @examples
#' data(example_seurat_obj)
#' example_seurat_obj
#' dim(example_seurat_obj)
#' head(example_seurat_obj@meta.data)
"example_seurat_obj"
