# Getting Started with Model Analysis

All the functions below assume you've already arranged your data in the format
expected by `mixOmics`.

## Model variance analysis

Interfold provides the `get.block.centroids()` function to extract the centroids of variance across blocks in a DIABLO model and export them as a data frame and plot.  The `get.model.variance()` function is also provided which will take either an sPLS-DA model or a DIABLO model as input and creates a data frame showing the percentage contributions of each component to the model variance.

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
