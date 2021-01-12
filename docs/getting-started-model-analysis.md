# Getting Started with Model Analysis

The following functions are provided to assist in the analysis and visualisation
of a model fitted by `mixOmics`.  Focus has been given to sPLS-DA models where a
single block is being fitted (single-omics) and to DIABLO models where multiple
blocks are being fitted simulataneously (multi-omics).  The code below will
generate known trained models and corresponding performance tests for each model
type, which are then used to demonstrate the functions in the rest of this
document.

```R
library(interfold.bio)
library(mixOmics)
library(dplyr)

#######################################
## sPLS-DA model and performance result
#######################################

# Prepare data
data(srbct)
splsda.X <- srbct$gene
splsda.Y <- srbct$class

# Model fit
splsda.ncomp <- 3
splsda.keepX <- c(9, 280, 30)
splsda.model <- splsda(splsda.X, splsda.Y,
                       ncomp = splsda.ncomp, keepX = splsda.keepX)

# Preformance result -- takes about 1 min to run
set.seed(1234) # for reproducibility
splsda.perf <- perf(splsda.model, validation = "Mfold", folds = 5,
                    dist = 'max.dist', nrepeat = 10,
                    progressBar = FALSE)


######################################
## DIABLO model and performance result
######################################

# Prepare data
data('breast.TCGA')
diablo.data <- list(mRNA = breast.TCGA$data.train$mrna, 
                    miRNA = breast.TCGA$data.train$mirna, 
                    proteomics = breast.TCGA$data.train$protein)
diablo.Y <- breast.TCGA$data.train$subtype
diablo.design <- matrix(0.1, ncol = length(diablo.data),
                        nrow = length(diablo.data), 
                        dimnames = list(names(diablo.data), names(diablo.data)))
diag(diablo.design) <- 0

# Model fit
diablo.ncomp <- 2
diablo.keepX <- list(mRNA = c(6, 14), miRNA = c(5, 18), proteomics = c(6, 7))
diablo.model <- block.splsda(X = diablo.data, Y = diablo.Y,
                             ncomp = diablo.ncomp, keepX = diablo.keepX,
                             design = diablo.design)

# Preformance result -- takes about 1 min to run
set.seed(4321) # for reproducibility
diablo.perf <- perf(diablo.model, validation = 'Mfold', M = 10, nrepeat = 10, 
                    dist = 'centroids.dist')
```

## Model variance analysis

InterFold provides the `get.block.centroids()` function to extract the centroids of variance across blocks in a DIABLO model and export them as a data frame and plot.  The `get.model.variance()` function is also provided which will take either an sPLS-DA model or a DIABLO model as input and creates a data frame showing the percentage contributions of each component to the model variance.

## Feature analysis for sPLS-DA models

- get.loadings.table()
- feature.selection.stability()
- merge.feature.stability()
- plot.feature.stability()

## Feature analysis for DIABLO models

- get.diablo.top.loadings()
- diablo.selection.stability()
- get.diablo.top.loadings.with.stability()
- get.diablo.top.features()
- find.feature.associations()
- export.matrix.as.network()

## Model predictivity

- plot.predicted.projection()

## Utility functions

- center.truncate()
