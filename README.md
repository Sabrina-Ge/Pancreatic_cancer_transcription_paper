This repository contains major scripts for the preprint:

**Extra-lineage tissue programs define the transcription states of human pancreatic cancer** 

Sabrina Ge, Paul Tonon, Jingxiong Xu, Gun Ho Jang, Ferris Nowlan, Jimin Min, Karen Ng, Eugenia Flores Figueroa, Amy Zhang, Michelle Chan-Seng-Yue, Yuanchang Fang, Adriana Migliorini, Julie M. Wilson, Anna Dodd, Sebastian Arcila-Barrera, Ayah Elqaderi, Ilinca Lungu, Yu Zhang, Stephanie Ramotar, Shawn Hutchinson, Daniela Bevacqua, Ayelet Borgida, Spring Holter, Pathum Kossinna, Ruth Isserlin, Veronique Voisin, Andreas Mund, Gary D. Bader, Maria Cristina Nostro, Erica S. Tsang, Robert C. Grant, Grainne M. O’Kane, David A. Tuveson, Oren Parnas, Federico Gaiti, Nina Steele, Jennifer J. Knox, Anirban Maitra, Hartland W. Jackson, Steven Gallinger, Faiyaz Notta

<https://www.biorxiv.org/content/10.64898/2026.04.29.721655v1>

# Repository Contents

- [ssp_bkr.R](./ssp_bkr.R): `predict.class.bkr()` - single sample classifier of the updated transcriptional classification scheme for bulk RNA sequencing
- [singler_subtype.R](./singler_subtype.R): `singler_subtype()` - correlation-based classifier for single cell RNA sequencing
- [penalized_module_score.R](./penalized_module_score.R): `add_penalized_module_score()` - scores query gene sets (modules) by building upon Seurat's `AddModuleScore` to favour broad expression of query genes
- [misc_functions.R](./misc_functions.R): functions for converting human to mouse homologues and generating program genelists
- [references](./references): reference files needed for running previous functions
- [scripts](/.scripts): code for recreating other major analyses in the paper
- [demo](/.demo): demo files needed to run `predict.class.bkr()`, `singler_subtype()`, and `add_penalized_module_score()` on a small subset

# Novel functions for classification and scoring

## Single sample classifier for bulk RNA-seq

[ssp_bkr.R](./ssp_bkr.R): `predict.class.bkr()`

`mat` mode outputs a list of gene pairs used in the classifier

```         
predict.class.bkr(mode = "mat")
```

`pred` (default) mode assigns a classification of `Classical1`, `Classical2`, `hybrid`, `Basal1`, or `Basal2` to each column (sample) of the gene x sample TPM matrix `mat`.

```         
predict.class.bkr(mat)
```

Other default parameters can be altered depending on the type of expression matrix input.


## Correlation-based classifier for scRNA-seq

[singler_subtype.R](./singler_subtype.R): `singler_subtype()`

Uses correlation-based SingleR algorithm (see: https://bioconductor.org/packages/release/bioc/html/SingleR.html) to score and assign classification to single cells using a reference derived from bulk data.

Inputs:

- `mat` Count matrix where genes (HGNC format) are rows and cells are columns. LogCPM values are recommended. Can be a sparse matrix.
- `reference_file` Path to reference bulk tumour RData containing the objects ref_mat (LogTPM HGNC gene by sample matrix) and ref_labels (named character vector, sample: subtype label). Default "bulk_reference.RData".
- `find_mixed` Classify cells as  `Mixed`. Default TRUE.

Outputs a list with the following:

- `$data` data.frame where rows are cells and columns represent the classification (`Classical1`, `Classical2`, `Mixed`, `Basal1`, or `Basal2`), normalized confidence score and normalized classification scores
- `$singler_pred` SingleR raw output

```
out <- singler_subtype(mat, reference_file="references/singler_reference.RData")
out$data
out$singler_pred

```

## Penalized gene set scoring for scRNA-seq

[penalized_module_score.R](./penalized_module_score.R): `add_penalized_module_score()`

Building on Seurat's `AddModuleScore`, calculates the average expression level of each query gene set per single cell when compared with randomly selected matched control features. The average expression is additionally penalized based on the fraction of genes whose expression exceeds the mean expression of matched control genes. Fewer genes with expression exceeding the control leads to greater penalization, allowing higher scores for more consistent query gene set expression.

Takes a Seurat object `dataset` and a named list of query genesets `genelists`, and outputs a Seurat object with scores added in the object meta data. 
The new metadata columns will be named after the names of the query geneset list.
Additional parameters change the penalty factor (`k`, `x0`) or the the number of expression bins for selecting the control genes (`nbin`) as in Seurat's `AddModuleScore`.

```
seurat_obj <- add_penalized_module_score(seurat_obj, genelists=named_list_of_genesets)
seurat_obj[[names(named_list_of_genesets)]]
```

# Demo

## Set-up

The demo analyses can be run on a standard computer with 16GB+ RAM and 4+ cores. We recommend using R version 4.3.0+. 
While the demo has been tested specifically on Ubuntu 20.04 and macOS Tahoe 26.5.2, R and all relevant packages should be compatible with most modern Windows, Mac, and Linux operating systems.

Install R and the development environment RStudio by following these instructions:

R: <https://cran.r-project.org/>

RStudio: <https://docs.posit.co/ide/user/>

Installation of R and Rstudio took around 2 minutes on an M5 MacBook Air with 24GB RAM running macOS Tahoe 26.5.2 at 1Gbps download.

## Dependencies

To install the dependencies, type the following code sections into an `R` session:

`predict.class.bkr()`: none

`singler_subtype()`: SingleR v2.4.1+ (from Bioconductor 3.18+)

```
install.packages("BiocManager")
BiocManager::install("SingleR") 
```

`add_penalized_module_score()`: Seurat v5.1.0+, dplyr v1.1.4+, ggplot2 v3.5.2+, Matrix v1.6-4+

```         
install.packages(c("dplyr", "Seurat", "ggplot2"))
install.packages("Matrix") # generally unneeded as Matrix usually comes pre-installed with R
```

Installation of dependencies took under 2 minutes on an M5 MacBook Air with 24GB RAM running macOS Tahoe 26.5.2 at 1Gbps download.


## Run and Results

Clone, or download and unzip the respository to your local computer. Open RStudio and using "File > Open Project...", open the `pancreatic_cancer_transcription_paper` directory as a project. In the "Files" tab, navigate to the `demo` directory and open `demo.R`. Source `demo.R` by selecting "Source" in the top right corner.

When successfully run, the script will produce output files related to each function in a new directory called `demo_outputs/`, which should match the existing files in `demo_output_example/`.

Running `demo.R` took under 1 minute on an M5 MacBook Air with 24GB RAM running macOS Tahoe 26.5.2.


