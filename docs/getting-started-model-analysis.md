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

InterFold provides the [`get.block.centroids()`](#get-block-centroids) function
to extract the centroids of variance across blocks in a DIABLO model and export
them as a data frame and plot.  The
[`get.model.variance()`](#get-model-variance) function is also provided which
will take either an sPLS-DA model or a DIABLO model as input and creates a data
frame showing the percentage contributions of each component to the model
variance.

### Get block centroids

Applying the `get.block.centroids()` function to the DIABLO model from above
gives the following:

```R
get.block.centroids(diablo.model)
## A tibble: 150 x 4
## Groups:   sample [150]
##    sample label comp1   comp2
##     <int> <fct> <dbl>   <dbl>
##  1      1 Basal  2.38  2.22  
##  2      2 Basal  1.74  2.03  
##  3      3 Basal  2.28 -0.0206
##  4      4 Basal  1.54  1.46  
##  5      5 Basal  1.93 -1.17  
## ... with 145 more rows
```

In addition to the tibble provided, a plot is produced showing the centroid
across the blocks in the first and second component for each sample in the data.
Each point in the plot represents a sample and is the mean of the three blocks
loaded for the DIABLO model.  The circles show the distribution of all the
points sharing the same label and gives a good indication of the separation of
classifications in each of the first two components.

![Block Centroids Plot](imgs/block-centroids-plot.png)

### Get model variance

The `get.model.variance()` function can be applied to the model to identify how
much variation each component captures in the model for a given block:

```R
get.model.variance(diablo.model, "mRNA")
##        Proportion Cumulative
## comp1 0.794465823  0.7944658
## comp2 0.008624793  0.8030906

get.model.variance(diablo.model, "miRNA")
##        Proportion Cumulative
## comp1 0.669987181  0.6699872
## comp2 0.001283394  0.6712706

get.model.variance(diablo.model, "proteomics")
##        Proportion Cumulative
## comp1 0.770439163  0.7704392
## comp2 0.002689493  0.7731287
```

Here we can see that the first component captures the vast majority of the
variance for all three blocks but that slightly less of the variance in `miRNA`
was captured than in `mRNA` or `proteomics`.

## Feature analysis for sPLS-DA models

InterFold provides functions for analysing the features selecting for fitting
sPLS-DA models.

- [`get.loadings.table()`](#get-feature-loadings-as-a-table) outputs a data
  frame of loadings for features selected in each component of the model.
- [`feature.selection.stability()`](#get-spls-da-feature-selection-stability)
  gives the loadings from above for a single component, alongside stability
  scores for the selected features.
- [`merge.feature.stability()`](#merge-feature-stability-on-a-table-of-loadings)
  takes a table such as that generated by `get.loadings.table()` and the output
  from `feature.selection.stability()` merging them into a single data frame.

### Get feature loadings as a table

By calling `get.loadings.table()` with the sPLS-DA model as a parameter, a data
frame will be produced listing the loadings of each feature selected for the
sparse feature set.  A separate column is included for each component and the
features are arranged according to their loading value in the earliest component
they are selected for.  Ordering is highest loading first.  Where a feature is
selected in more than one component, it appears only once in the table, showing
both loadings values.  Where features are not selected for a component, the
loading value is `NA`.

The code example below demonstrates typical output from this function.  The
output has been truncated for brevity in this getting started guide.  The actual
output is 320 lines long.

```R
get.loadings.table(splsda.model)
##            comp1         comp2       comp3
## g123  0.64922048            NA          NA
## g846  0.44828969            NA          NA
## g1606 0.30641602            NA          NA
## ...
## g1389         NA -2.312917e-01          NA
## g246          NA -2.001212e-01          NA
## g545          NA -1.986890e-01          NA
## ...
## g340          NA -4.226175e-02          NA
## g742          NA  4.204118e-02  0.14003503
## g1143         NA -4.124087e-02          NA
## ...
## g174          NA            NA -0.35418530
## g1896         NA            NA -0.29939099
## g603          NA            NA -0.27261886
```

### Get sPLS-DA feature selection stability

Calling `feature.selection.stability()` with the model, performance test result
and a component number as parameters generates a data frame showing the loadings
values are were returned above, but for a single component only.  The second
column in the data frame gives the stability value of the feature.  This is
expressed as a fraction and indicates how many of the models in the performance
test included this feature in the sparse set.  High stability indicates that a
feature is very important to the fit of the model and is therefore probably
important to separating the classes the model is trying to discriminate.

The example below shows extraction of the stability values for the 9 top
features included in component 1.  Note how all stabilities are an exact
multiple of 0.02.  This is because our performance test used 10 repeats and 5
folds, so 50 models were fitted, and each time a model selects a feature, this
adds 0.02 to the stability of that feature.

```R
feature.selection.stability(splsda.model, splsda.perf, 1)
##            value stability
## g123  0.64922048      0.50
## g846  0.44828969      0.50
## g1606 0.30641602      0.40
## g335  0.30031495      0.48
## g836  0.26392532      0.46
## g783  0.22022625      0.28
## g758  0.22019689      0.38
## g1386 0.16248889      0.38
## g585  0.02058632      0.12
```

### Merge feature stability on a table of loadings

After getting the loadings table for all the components of the model, it is
probably desirable to create a table with stability values shown in a final
column.  `merge.feature.stability()` can be called with the loadings table and
stability table as parameters.  The returned data frame adds a column for the
stability values matching the features in the loadings table.  `NA` is used
where the stability value for a feature was not in the stability table.

Output for the loadings table is truncated in the example below for brevity.
You will see the full output in your own R instance.

```R
# Prepare the loadings table and stability value table we intend to merge on.
loadings.table <- get.loadings.table(splsda.model)
stability.comp3 <- feature.selection.stability(splsda.model, splsda.perf, 3)

# Check the raw loadings table before merging any stability values.
loadings.table
##            comp1         comp2       comp3
## g123  0.64922048            NA          NA
## g846  0.44828969            NA          NA
## g742          NA  4.204118e-02  0.14003503
## g1143         NA -4.124087e-02          NA
## g174          NA            NA -0.35418530
## g1896         NA            NA -0.29939099

# Merge on stability values for component 3.
merge.feature.stability(loadings.table, stability.comp3)
## feature      comp1         comp2       comp3 stability
##    g123 0.64922048            NA          NA        NA
##    g846 0.44828969            NA          NA        NA
##    g742         NA  4.204118e-02  0.14003503      0.76
##   g1143         NA -4.124087e-02          NA        NA
##    g174         NA            NA -0.35418530      0.34
##   g1896         NA            NA -0.29939099      0.32
```

## Feature analysis for DIABLO models

InterFold provides a number of functions for analysing the output from DIABLO
models produced in `mixOmics`.

- [`get.diablo.top.loadings()`](#get-diablo-top-loadings) outputs a data frame
  of feature loadings for a specified block in the DIABLO model.
- [`diablo.selection.stability()`](#get-diablo-feature-selection-stability)
  extracts a data frame of stability scores from a DIABLO model performance test
  for features in a particular block and component.
- [`get.diablo.top.loadings.with.stability()`](#get-diablo-top-loadings-with-stability)
  provides a concatenated data frame of feature loadings and stability scores
  for a given block in the DIABLO model.
- [`get.diablo.top.features()`](#get-diablo-top-features) applies a correction
  factor to stability scores so that features can be ranked across different
  blocks and different components, returning a sorted data frame of features in
  the model.
- [`plot.feature.stability()`](#plot-feature-stability) provides a simple
  histogram plot of the highest stability features in a single component and
  block, allowing visual comparison of their stability scores.
- [`find.feature.associations()`](#find-diablo-feature-associations)
- [`export.matrix.as.network()`](#export-diablo-matrix-as-a-network)

### Get DIABLO top loadings

`get.diablo.top.loadings()`

```R
get.diablo.top.loadings(diablo.model$loadings$proteomics)
##                   comp1       comp2
## ER-alpha    -0.74995295          NA
## GATA3       -0.62448770          NA
## ASNS         0.16825721          NA
## Cyclin_B1    0.12026614          NA
## AR          -0.06842432 -0.11183919
## JNK2        -0.01137347          NA
## HER2                 NA -0.67182918
## HER2_pY1248          NA -0.62130938
## EGFR_pY1068          NA -0.33229248
## c-Kit                NA  0.17791016
## HER3_pY1289          NA -0.08706226
## XRCC1                NA  0.02149536
```

### Get DIABLO feature selection stability

`diablo.selection.stability()`

```R
diablo.selection.stability(diablo.perf, comp = 1, block = 'proteomics')
##             feature stability
## ER-alpha   ER-alpha      1.00
## GATA3         GATA3      1.00
## ASNS           ASNS      0.99
## Cyclin_B1 Cyclin_B1      0.98
## AR               AR      0.80
## JNK2           JNK2      0.54
## PR               PR      0.41
## Cyclin_E1 Cyclin_E1      0.22
## INPP4B       INPP4B      0.12
```

### Get DIABLO top loadings with stability

`get.diablo.top.loadings.with.stability()`

```R
get.diablo.top.loadings.with.stability(diablo.model, diablo.perf, block = 'proteomics', feature.count = 50)
##                   comp1       comp2 stability
## ER-alpha    -0.74995295          NA      1.00
## GATA3       -0.62448770          NA      1.00
## ASNS         0.16825721          NA      0.99
## Cyclin_B1    0.12026614          NA      0.98
## AR          -0.06842432 -0.11183919      0.80
## JNK2        -0.01137347          NA      0.54
## HER2                 NA -0.67182918      1.00
## HER2_pY1248          NA -0.62130938      1.00
## EGFR_pY1068          NA -0.33229248      1.00
## c-Kit                NA  0.17791016      1.00
## HER3_pY1289          NA -0.08706226      0.95
## XRCC1                NA  0.02149536      0.56
```

### Get DIABLO top features

`get.diablo.top.features()`

```R
get.diablo.top.features(diablo.model, diablo.perf, feature.count = 150)
## rank rank.score      block component         feature stability
##    1 1.00000000      miRNA         1    hsa-mir-130b 1.0000000
##    2 1.00000000      miRNA         1      hsa-mir-17 1.0000000
##    3 1.00000000      miRNA         1     hsa-mir-505 1.0000000
##    4 1.00000000      miRNA         1     hsa-mir-590 1.0000000
##    5 1.00000000       mRNA         1           KDM4B 1.0000000
##    6 1.00000000       mRNA         1          ZNF552 1.0000000
##    7 1.00000000 proteomics         1        ER-alpha 1.0000000
##    8 1.00000000 proteomics         1           GATA3 1.0000000
##    9 1.00000000 proteomics         2           c-Kit 1.0000000
##   10 1.00000000 proteomics         2     EGFR_pY1068 1.0000000
## ...
##  109 0.05698361      miRNA         2     hsa-mir-23b 0.1000000
##  110 0.05698361      miRNA         2    hsa-let-7a-1 0.1000000
##  111 0.05698361      miRNA         2     hsa-mir-30d 0.1000000
##  112 0.05698361      miRNA         2     hsa-mir-451 0.1000000
##  113 0.05698361      miRNA         2     hsa-mir-144 0.1000000
```

### Plot feature stability

`plot.feature.stability()`

```R
stability <- diablo.selection.stability(diablo.perf, comp = 2, block = 'miRNA')
plot.feature.stability(stability)
```

![Feature Stability Plot](img/../imgs/plot-feature-stability.png)

### Find DIABLO feature associations

`find.feature.associations()`

```R
find.feature.associations(diablo.model, block.count = 3)
##                   KDM4B     ZNF552       FUT8      LRIG1      CCNA2      PREX1
## KDM4B         1.0000000  0.8247985  0.5777819  0.6780996 -0.6194948  0.7235040
## ZNF552        0.8247985  1.0000000  0.6428172  0.7072751 -0.6399958  0.7308222
## FUT8          0.5777819  0.6428172  1.0000000  0.5162760 -0.4628657  0.5168250
## LRIG1         0.6780996  0.7072751  0.5162760  1.0000000 -0.5223886  0.5995716
## CCNA2        -0.6194948 -0.6399958 -0.4628657 -0.5223886  1.0000000 -0.5468418
## PREX1         0.7235040  0.7308222  0.5168250  0.5995716 -0.5468418  1.0000000
```

![Feature Associations Plot](imgs/find-feature-associations.png)

### Export DIABLO matrix as a network

`export.matrix.as.network()`

```R
associations <- find.feature.associations(diablo.model, block.count = 3)
export.matrix.as.network(associations, filename = "network.csv", cutoff = 0.7,
                         block.feature.count = 17)
##    feature.1    feature.2      value block
##        KDM4B       ZNF552  0.8247985     1
##        KDM4B        PREX1  0.7235040     1
##       ZNF552        LRIG1  0.7072751     1
##       ZNF552        PREX1  0.7308222     1
```

## Model predictivity

`plot.predicted.projection()`

```R
# Prepare data
data(srbct)
all.X <- srbct$gene
all.Y <- srbct$class
all.index <- 1:length(all.Y)

# Subset the data for a training set
set.seed(1234) # for reproducibility of subsetting
train.index <- sort(sample(all.index, length(all.index) * 0.8))
train.X <- all.X[train.index,]
train.Y <- all.Y[train.index]

# Model fit with the training set
predict.model <- splsda(train.X, train.Y, ncomp = 3, keepX = c(9, 280, 30))

# Prepare the prediction set
predict.index <- setdiff(all.index, train.index)
predict.X <- all.X[predict.index,]
predict.Y <- all.Y[predict.index]

# Create a prediction for the predict set and plot using plot.predicted.projection()
prediction <- predict(predict.model, predict.X)
projection.plot <- plot.predicted.projection(prediction, predict.Y)
print(projection.plot)
```

![Prediction Projection Plot](imgs/plot-predicted-projection.png)

## Utility functions

`center.truncate()`

```R
center.truncate("When a string is particularly long it will be truncated back to 43 characters")
## "When a string is par...back to 43 characters"

labels = c("Short name",
           "Very long name would need truncating for plotting",
           "Plots don't work well when the axis is forced over by long labels")
unname(sapply(labels, center.truncate))
## [1] "Short name"
## [2] "Very long name would...uncating for plotting"
## [3] "Plots don't work wel...d over by long labels"
```
