# Omicsfold

![Maturity level-Prototype](https://img.shields.io/badge/Maturity%20Level-Prototype-red)

![](omicsfold_id.png)

### Multi-omics data normalisation, model fitting, and visualisation.

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

See the [`DESCRIPTION`](OmicsFold/DESCRIPTION) file for a complete
dependency list

### Installation

Due to the number of dependencies and the number of places those dependencies
come from, there is an installation script available.  This can be run by
opening up an R session in your preferred environment, ensuring your working
directory is the `OmicsFold` directory, then issuing the following commands:

```R
source('install.R')
install.omicsfold()
```

This should install all the dependencies and then finally the OmicsFold package
itself.  If there are any issues due to versions changing or changes in which
repository maintains the active version of a package, you may have to update the
script.

If you are having issues installing OmicsFold in a conda environment, please try
the following steps: 

First, create the conda environment:
```Shell
conda create --name OmicsFold 
source activate OmicsFold
conda install r=3.6.0
conda install -c conda-forge boost-cpp
```

Second, launch R in the conda environment and manually install the following packages (or if you are installing directly in a local instance of R):
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("metagenomeSeq")
BiocManager::install("org.Mm.eg.db")
install.packages("XML", repos = "http://www.omegahat.net/R")
source("http://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/GeneAnnotation/installAnRichment.R")
installAnRichment()
source('install.R')
install.omicsfold()
```
For installation using nextflow (https://www.nextflow.io/docs/latest/getstarted.html) please see https://github.com/AstraZeneca/Omicsfold/tree/master/OmicsFold/nextflow_pipeline

### Usage

Import the `OmicsFold` and the `mixOmics` packages in R and you're ready to
go.  Some functions also require `dplyr` to be loaded so it's a good idea to
load it anyway. Certain plotting functions also may require ggplot2 to be loaded.

```R
library(OmicsFold)
library(mixOmics)
library(dplyr)
library (ggplot2) #(optional)
```

### Data Normalisation

A number of normalisation functions have been provided.  Each has documentation
which can be read in the usual way in R.  For example, the help for the function
`normalise.tss` can be viewed by calling `?normalise.tss`.  A brief description
of the usage of each function can be read in the [Getting Started with
Normalisation](docs/getting-started-normalisation.md) document, with a few key
functions also showing example code for how to use it.

- `low.count.removal()`
- `normalise.tss()`
- `normalise.css()`
- `normalise.logit()`
- `normalise.logit.empirical()`
- `normalise.clr()`
- `normalise.clr.within.features()`

### Analysis of mixOmics Output

Once a `mixOmics` model has been fitted, OmicsFold can be used to perform a
number of visualisation and data extraction functions.  Below is a brief list of
the functionality provided.  While these are well documented in the R help
system, descriptions of how to use each function can also be found in the
[Getting Started with Model Analysis](docs/getting-started-model-analysis.md)
document.

- **Model variance analysis** - functions are provided to extract the percentage
  contributions of each component to the model variance and the centroids of
  variance across the blocks of a DIABLO model.
- **Feature analysis for sPLS-DA models** - feature loadings on the fitted
  singleomics model can be exported as a sorted table, while feature stability
  across many sparse model fits can also be exported.  As there may be many
  components to export stability for, another function lets you combine these
  into a single table as well as a plotting function allowing you to plot
  stability of the selected features as a visualisation.
- **Feature analysis for DIABLO models** - similarly to the features for
  singleomics models above, multiomics models can also have feature loadings and
  stability exported. Associated correlations between features of different 
  blocks can be exported as either a matrix and then also converted to a CSV 
  file appropriate for importing into Cytoscape where it can form a network 
  graph.
- **Model predictivity** - we provide a function to plot the predictivity of a
  model from a confusion matrix.
- **Utility functions** - offers a way to take long feature names being passed
  to plots and truncate them for display.
- **BlockRank** - implements a novel approach to analysing feature importance 
  between blocks of data.



## Other Information

To contact the maintainers or project director, please refer to the
[`AUTHORS`](AUTHORS.md) file.  If you are thinking of contributing to OmicsFold,
all the information you will need is in the [`CONTRIBUTING`](CONTRIBUTING.md)
file.

OmicsFold is licensed under the [Apache-2.0 software
licence](https://www.apache.org/licenses/LICENSE-2.0) as documented in the
[`LICENCE`](LICENCE.md) file.  Separately installed dependencies of OmicsFold
may be licensed under different licence agreements.  If you plan to create
derivative works from OmicsFold or use OmicsFold for commercial or profitable
enterprises, please ensure you adhere to all the expectations of these
dependencies and seek legal advice if you are unsure.
