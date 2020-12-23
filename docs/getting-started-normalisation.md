# Getting Started with Normalisation

All the functions below assume you've already arranged your data in the format
expected by `mixOmics`.  In other words, you should have samples in rows and
features in columns, as either a data frame, matrix or other table like data
structure.

Data loaded to demonstrate each function is demo data from the `mixOmics`
package so that should be installed and loaded in order to try out these
functions.

## Low Count Removal

`low.count.removal()` removes features from the data which are unlikely to contribute to the fit of a
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

## Total Sum Scaling

`normalise.tss()` normalises count data sample-by-sample, to a scale of 0..1, using Total Sum
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

## Cumulative Sum Scaling

## Logit

There is also a variant for applying the logit transformation in an empirical way.

## Centered Log Ratio

There is also a variant of the centered log ratio for within features normalisation.
