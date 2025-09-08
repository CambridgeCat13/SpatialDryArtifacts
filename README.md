
# SpatialDryArtifacts

<!-- badges: start -->
<!-- badges: end -->

## Overview

SpatialDryArtifacts implements enhanced methods for detecting edge dryspots - technical artifacts caused by incomplete reagent coverage during Visium spatial transcriptomics sample preparation. This package extends existing spatial transcriptomics quality control methods by incorporating spatial neighborhood information.

## Key Features

- Detection of edge dryspots using raster-based focal transformations
- Classification of problem areas into removal categories
- Integration with existing Bioconductor workflows
- Compatibility with SpatialExperiment objects
  
## Installation

You can install the development version of SpatialDryArtifacts from [GitHub](https://github.com/) with:

```r
# install.packages("pak")
pak::pak("CambridgeCat13/SpatialDryArtifacts")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SpatialDryArtifacts)

# Detect edge dryspots
spe <- detectEdgeDryspots(spe, qc_metric = "sum_gene")

# Classify results  
spe <- classifyEdgeDryspots(spe)
```

