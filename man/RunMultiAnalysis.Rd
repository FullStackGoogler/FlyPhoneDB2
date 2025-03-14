% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunFlyPhone.R
\name{RunMultiAnalysis}
\alias{RunMultiAnalysis}
\title{Multi Analysis}
\usage{
RunMultiAnalysis(
  knowledgebase_version,
  counts_fn = NULL,
  metadata_fn = NULL,
  delimitor = "\\t",
  seuratObject = NULL,
  DEG = NULL,
  control_name = NULL,
  mutant_name = NULL,
  pct_filter = 0.1,
  perm_times = 1000,
  deletePE = TRUE,
  base_output_dir = ""
)
}
\arguments{
\item{knowledgebase_version}{Version of knowledgebase to use. Valid values are: "Version 1", "Version 2 All", "Version 2 High", and "Version 2 High/Moderate".}

\item{counts_fn}{The filepath to the counts matrix.}

\item{metadata_fn}{The filepath to the accompanying metadata file.}

\item{delimitor}{The separator for the counts file if it is a .txt file. Default is tab.}

\item{seuratObject}{Optionally can upload a Seurat Object instead of a matrix/metadata file pair.}

\item{DEG}{The filepath to the DEG file, if provided.}

\item{control_name}{The name of the control sample, if applicable.}

\item{mutant_name}{The name of the mutant sample, if applicable.}

\item{pct_filter}{Percent threshold to filter expression of genes by. Default value of 0.1.}

\item{perm_times}{The number of times to shuffle cluster alignments and calculate ligand/receptor averages for p-value calculation in raw interaction scoring. Default value of 1000.}

\item{deletePE}{Whether or not to delete the percentage expression file after the pipeline is finished. Defaults to true.}

\item{base_output_dir}{The directory for FlyPhone to send all result files to. Defaults to the current working directory.}
}
\description{
Generate multi sample analysis data.
}
