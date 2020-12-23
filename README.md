# InterFold

![Maturity level-0](https://img.shields.io/badge/Maturity%20Level-ML--0-red)

Multi-omics data normalisation, model fitting and visualization.

## Overview

This is a utility R package containing custom code and scripts developed to
establish a working approach for integration of multi-omics data.

The package provides a unified toolkit for the analysis and integration of
multi-omic high-throughput data. It relies upon the
[`mixOmics`](http://mixomics.org/) toolkit to provide implementations of many of
the underlying projection to latent structures (PLS) methods used to analyse
high-dimensional data. In addition to this, it includes custom implementations
of data pre-processing, normalisation, collation, model validation,
visualisation & output functions.

The originally individual scripts have been collected into a formal package that
should be installable and usable within an analysts' R environment without
further configuration. The package is fully documented at the function level.

## Getting Started

This package and analysis requires R v3.6 or above. It is largely built upon the
`mixOmics` integration framework. The dependencies vary significantly in source,
so an installation script is provided to make satisfying the dependencies as
simple as possible. `mixOmics` installs its own dependencies as well. Note that
we install `mixOmics` from the GitHub repository as this version is more up to
date than the one on Bioconductor and has a number of fixes which are needed to
avoid bugs.

Notable dependencies that will be installed if they are not already:

- mixOmics
- WGCNA
- ggplot2
- dplyr & magrittr
- reshape2

See `DESCRIPTION` for a complete dependency list

### Installation

Due to the number of dependencies and the number of places those dependencies
come from, there is an installation script available.  This can be run by
opening up an R session in your preferred environment, ensuring your working
directory is the `interfold.bio` directory, then issuing the following commands:

```R
source('install.R')
install.interfold.bio()
```

This should install all the dependencies and then finally the InterFold package
itself.  If there are any issues due to versions changing or changes in which
repository maintains the active version of a package, you may have to update the
script.

### Usage

Import the InterFold and the mixOmics packaege and you're ready to go.

```R
library(interfold.bio)
library(mixOmics)
```

### Data Normalisation

### Visualization of mixOmics Output

### WGCNA Analysis of Expression Data
