
<!-- README.md is generated from README.Rmd. Please edit that file -->

# omics4drug <img src="man/figures/logo.png" width = "175" height = "200" align="right" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/yen-kim/omics4drug/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yen-kim/omics4drug/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![DOI](https://zenodo.org/badge/1056782986.svg)](https://doi.org/10.5281/zenodo.17117623)
<!-- badges: end -->

An R toolkit to facilitate Mass Spectrometry-based Proteomics and
Phosphoproteomics data analysis

`omics4drug` is designed for the analysis and visualization of Mass
Spectrometry-based phosphoproteomics and proteomics data in drug
discovery. The package provides functions for quality control,
normalization, pathway enrichment analysis, and drug-target prediction.

## Installation

To get the latest in-development features, install the development
version from GitHub:

``` r
if(!requireNamespace("devtools", quietly = TRUE)) {
 install.packages("devtools")
}
devtools::install_github("yen-kim/omics4drug")
```

This package is also accessible for download via Zenodo with the DOI
[10.5281/zenodo.17117624](https://doi.org/10.5281/zenodo.17117623).

### Functions

See [Package
index](https://yen-kim.github.io/omics4drug/reference/index.html) for
full list of functions.

1.  Data Processing and Quality Control

- `get_count_phosphosite()`: Counts and visualizes the number of unique
  phosphosites per sample or group, often based on a probability
  threshold.  
- `get_count_protein()`: Counts and visualizes the number of unique
  protein groups per sample or group.  
- `get_cv()`: Calculates and visualizes the coefficient of
  variation (CV) for a given dataset, useful for assessing data
  variability and quality.
- `get_sty()`: Calculates and visualizes the count and percentage of
  phosphorylation sites (Serine (S), Threonine (T), Tyrosine (Y)).

2.  Data Normalization

- `get_norm_phos()`: Normalizes phosphosite intensity data to account
  for variations between samples.
- `get_norm_prot()`: Normalizes protein group intensity data.

3.  Functional and Pathway Enrichment Analysis

- `get_GO()`: Performs Gene Ontology (GO) enrichment analysis to
  identify biological processes, molecular functions, or cellular
  components that are overrepresented in your data.
- `get_KEGG()`: Performs KEGG pathway enrichment analysis to determine
  which biological pathways are significantly impacted.

4.  Kinase and Drug Prediction

- `get_KSEA()`: Performs Kinase Substrate Enrichment Analysis (KSEA) to
  predict the activity of kinases based on the phosphorylation of their
  substrates.
- `get_inhibitor()`: Predicts which drugs might target the kinases
  identified in your analysis, using an external database.

5.  Others

- `get_annotation()`: Map Gene Identifiers

## Acknowledgements

This R package was produced with support from the [Copenhagen
University](https://www.ku.dk/en) through the [DISCOVER PhD
program](https://discover.ku.dk).
