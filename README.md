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

Import the InterFold and the mixOmics package and you're ready to go.

```R
library(interfold.bio)
library(mixOmics)
```

### Data Normalisation

A number of normalisation functions have been provided.  Each has documentation
which can be read in the usual way in R.  For example, the help for the function
`normalise.tss` can be viewed by calling `?normalise.tss`.  A brief description
of the usage of each function made available is below, with a few key functions
also showing a demo using data supplied with the `mixOmics` package.

All the functions below assume you've already arranged your data in the format
expected by `mixOmics`.  In other words, you should have samples in rows and
features in columns, as either a data frame, matrix or other table like data
structure.

#### `low.count.removal`

Removes features from the data which are unlikely to contribute to the fit of a
model because they show low counts/expression relative to the rest of the data.
The higher the percentage provided, the more features will be discarded.

```R
data(Koren.16S)
dim(Koren.16S$data.raw)
## [1]  43 980

normalised <- low.count.removal(Koren.16S$data.raw, 0.03)
dim(normalised)
## [1]  43 816
```

#### `normalise.tss`

Normalises count data sample-by-sample, to a scale of 0..1, using Total Sum
Scaling.  This accounts for sequencing differences between samples.  After this
transformation, all samples will sum to 1.0 and values for each feature will be
relative.  Values can be offset from zero by providing the optional `offset`
argument.  In the example below, we compare the TSS function from InterFold with
the pre-normalised TSS data in the Koren 16S data set.

```R
data(Koren.16S)
Koren.16S$data.TSS[1:3, 3:5]
##                410908       177792      4294607
## Feces659 0.0002961208 0.0293159609 0.0002961208
## Feces309 0.0003447087 0.0003447087 0.0003447087
## Mouth599 0.0004083299 0.0002041650 0.0002041650

# Now apply our own TSS normalisation to the raw data
normalised <- normalise.tss(Koren.16S$data.raw)
normalised[1:3, 3:5]
##                410908       177792      4294607
## Feces659 0.0002961208 0.0293159609 0.0002961208
## Feces309 0.0003447087 0.0003447087 0.0003447087
## Mouth599 0.0004083299 0.0002041650 0.0002041650
```

#### `normalise.css`

#### `normalise.logit.empirical`

#### `normalise.logit`

#### `normalise.clr`

#### `normalise.clr.within.features`

### Visualization of mixOmics Output

### WGCNA Analysis of Expression Data
