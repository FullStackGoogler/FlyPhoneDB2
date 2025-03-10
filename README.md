# FlyPhoneDB2

## Overview

FlyPhone is a tool for analyzing cell-cell communications in *Drosophila Melanogaster* (fruit fly) single-cell RNA-sequencing datasets.

## Changes in V2

- New functionality to investigate the differences of cell-cell communication events for multi-samples.
- Core algorithm has been improved, resulting in significantly faster compute time.
- Ligand-Receptor database has been greatly expanded to include 1,774 pairs with varying degrees of confidence, ranging between "Low", "Moderate", and "High".
- Added more cell-cell communication visualizations, as well as improving on previous visualization types.
- Added visualizations to analyze pathway activity by taking into consideration downstream reporter genes.

## Installation
Installation is simple with the use of the [devtools](https://devtools.r-lib.org/) package.

```R
# devtools package should be installed and loaded prior to running this command
install_github('FullStackGoogler/FlyPhoneDB2')
```

## Usage

The entire FlyPhone pipeline can be ran with `RunFlyPhone()`. See below for detailed documentation on each of its parameters.

```R
RunFlyPhone(
    knowledgebase_version,
    counts_fn,
    metadata_fn,
    delimitor,
    seuratObject,
    DEG,
    control_name,
    mutant_name,
    pct_filter,
    perm_times,
    deletePE,
    base_output_dir
)
```

### Input Data

There are two methods for uploading your scRNA-seq data:

1. Upload a counts matrix and metadata file
    - If there are multiple samples, there should be a column labeled `Condition`.
    - Metadata should always be a `.CSV` file. Counts can either be a `.CSV` file or a `.TXT` file; if it is a `.TXT` file, the `delimitor` argument should be specified.

2. Upload a Seurat Object **(Recommended)** 
    - Recommended data format, as large Counts data files can take a long time to be read in and processed.
    - The meta.data in the Seurat Object should still have a column labeled `Condition`.

3. Differentially Expressed Genes (DEG) file
    - This file is required for the multi-sample DEG analysis. If there are multiple samples but a DEG file is not provided, FlyPhone will run the equivalent of single-sample analysis on each sample. The required columns are listed below:
        1. `Gene` - Gene name
        2. `avg_log2FC` - Log 2 Fold Change value
        3. `cell_type` - The cell type associated with a given gene

### Required Parameters

These variables are required, regardless of the data uploaded:
- `knowledgebase_version` - Specifies which version of FlyPhone's knowledge base (Ligand-Receptor pairs & Core components information) to use when running the pipeline. Valid options are:

    1. `Version 1` - Utilizes FlyPhone v1's version
    2. `Version 2 All` - Utilizies FlyPhone v2's version
    3. `Version 2 High` - Filters the L-R pairs to include those with a "High" rank
    4. `Version 2 High/Moderate` - Filters the L-R pairs to include those with either a "High" or "Moderate" rank

These variables are required if there are multiple samples present
- `control_name` - The control condition name; value should match that in the `Condition` column in your metadata.
- `mutant_name` - The mutant condition name; value should match that in the `Condition` column in your metadata.

### Optional Parameters

The remaining arguments are optional and have default values:

- `pct_filter` - Any genes with a percent expression across all cells below this threshold will be ignored in the analysis. Default value is 10% (0.1).
- `perm_times` - The number of times to run the permutation test for generating the p-value scores that determine the significance of a Ligand-Receptor signaling pair for any given celltypes. Larger values may result in more accurate results, but will take longer to run. Default value is 1,000.
- `deletePE` - `RunFlyPhone()` generates a percentage expression file for later analysis. If desired, set this argument to `FALSE` to keep the file (which can be found in the `.temp` folder); otherwise, it will be deleted upon completion. Default value is `TRUE`.
- `base_output_dir` - The output directory path for of `RunFlyPhone()`'s. Default value is `""`, or the current working directory.

## Web-based version

In addition to the R package, we've also built a website with the following functionality:
- **Analysis Pipeline -** GUI for the FlyPhone pipeline. Users can upload their data, where it will then be ran on the cloud. Once finished, you'll have a link containing your results.
- **Annotation Explorer -** Explore the latest version of our knowledge base (currently v2.5).
- **Community Annotation -** A form that allows you to submit your own Ligand-Receptor annotation, where it will then be reviewed before being added to the knowledge base.
- **Example Analysis -** To provide users with a better idea of what our tool returns, we've analyzed a fullbody scRNA-seq dataset from Yorkie tumor flies and have provided explanations.

**Link:** https://flyrnai.org/tools/fly_phone_2/web/

![image](https://github.com/user-attachments/assets/f71578ec-91a6-405e-a565-656a0a045dbd)

## Help/Contact

//TODO: Finish
