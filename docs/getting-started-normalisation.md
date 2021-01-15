# Getting Started with Normalisation

All the functions below assume you've already arranged your data in the format
expected by `mixOmics`.  In other words, you should have samples in rows and
features in columns, as either a data frame, matrix or other table like data
structure.

Data loaded to demonstrate each function is demo data from the `mixOmics`
package so that should be installed and loaded in order to try out these
functions.

## Low Count Removal

`low.count.removal()` removes features from the data which are unlikely to
contribute to the fit of a model because they show low counts/expression
relative to the rest of the data. The higher the percentage provided, the more
features will be discarded.

```R
data(Koren.16S)
dim(Koren.16S$data.raw)
## [1]  43 980

normalised <- low.count.removal(Koren.16S$data.raw, 0.03)
dim(normalised)
## [1]  43 816
```

## Total Sum Scaling

`normalise.tss()` normalises count data sample-by-sample, to a scale of 0..1,
using Total Sum Scaling.  This accounts for sequencing differences between
samples.  After this transformation, all samples will sum to 1.0 and values for
each feature will be relative.  Values can be offset from zero by providing the
optional `offset` argument.  In the example below, we compare the TSS function
from OmicsFold with the pre-normalised TSS data in the Koren 16S data set.

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

`normalise.css()` applies cumulative sum scaling normalisation to count data for
inter-sample depth.  This is an alternative to using total sum scaling and
relies on the implementation by `metagenomeSeq`.  It reformats `mixOmics` input
data so that it can be processed by `metagenomeSeq` and then converts the CSS
normalised output back to `mixOmics` input.  The definition for this
normalisation approach according to `metagenomeSeq` is as follows:

> Calculates each column's quantile and calculates the sum up to and including
> that quantile.

Below is an example of applying this normalisation to the same Koren 16S data
set as was used in the TSS example above:

```R
data(Koren.16S)
Koren.16S$data.raw[1:3, 3:5]
##            410908   177792  4294607
## Feces659        1       99        1
## Feces309        1        1        1
## Mouth599        2        1        1

# Now apply our CSS normalisation to the raw data
normalised <- normalise.css(Koren.16S$data.raw)
normalised[1:3, 3:5]
##            410908   177792  4294607
## Feces659 1.187222 6.993638 1.187222
## Feces309 1.179016 1.179016 1.179016
## Mouth599 1.711633 1.096030 1.096030
```

Here we see that the lowest counts of 1 for each feature / sample have much less
variance under CSS scaling, when compared to TSS scaling.

## Logit

`normalise.logit()` provides normalisation based on the [logit
function](https://en.wikipedia.org/wiki/Logit) which transforms 0.5 to zero,
values below 0.5 become negative and above 0.5 become positive.  The scale of
that negative or positive value is exponential and reaches negative/positive
infinity at 0.0 and 1.0 respectively.  This can be a useful transformation for
values in the 0..1 scale, bringing them back into Euclidean space after TSS
normalisation.  Below is an example of transforming values in this way.

```R
data(Koren.16S)
Koren.16S$data.TSS[1:3, 3:5]
##                410908       177792      4294607
## Feces659 0.0002961208 0.0293159609 0.0002961208
## Feces309 0.0003447087 0.0003447087 0.0003447087
## Mouth599 0.0004083299 0.0002041650 0.0002041650

# Now apply our logit normalisation to the TSS data
normalised <- normalise.logit(Koren.16S$data.TSS)
normalised[1:3, 3:5]
##                410908       177792      4294607
## Feces659    -8.124447    -3.499869    -8.124447
## Feces309    -7.972466    -7.972466    -7.972466
## Mouth599    -7.803027    -8.496378    -8.496378
```

As can be seen, this adds more distance between values, which can be more
beneficial for the model fitting.  In addition, values which are virtually zero
will be heavily modified towards a very negative value.  If any values to be
transformed are actually at 0.0 or 1.0 the logit function will generate infinity
values, which are inappropriate for modelling.  For this reason, a second
empirical function is provided by OmicsFold, `normalise.logit.empirical()`,
which moves measurements away from 0.0 and 1.0 on a per-feature basis, avoiding
the generation of infinity values.

## Centered Log-Ratio

`normalise.clr()` applies the centered log-ratio (CLR) transformation to the
data, where each measurement is divided by the mean of all the measurements and
the log of that ratio returned.  The ratio centers the mean of the data about
zero after the log is applied.  This transformation can be used with sum scaled
OTU data as an alternative to the logit transformation above.

This CLR function is a convenience wrapper around the
[`logratio.transfo()`](https://www.rdocumentation.org/packages/mixOmics/versions/6.3.2/topics/logratio.transfo)
function from the `mixOmics` package, maintaining the structure of the input
data after the transformation is applied.  Logs of zero produce infinity values,
so if the data being normalised contains zero values, a small offset can be
provided to prevent this problem.  The following is an example of applying the
CLR transformation with and without a small offset to the Koren 16S TSS data.
The offset is unnecessary in this data as there are no zero values, but the
results are similar for non-zero values.

```R
data(Koren.16S)
Koren.16S$data.TSS[1:3, 3:5]
##                410908       177792      4294607
## Feces659 0.0002961208 0.0293159609 0.0002961208
## Feces309 0.0003447087 0.0003447087 0.0003447087
## Mouth599 0.0004083299 0.0002041650 0.0002041650

# Now apply our CLR normalisation to the TSS data
normalised.1 <- normalise.clr(Koren.16S$data.TSS)
normalised.1[1:3, 3:5]
##                410908       177792      4294607
## Feces659   -0.3861923    4.2089276   -0.3861923
## Feces309   -0.3387886   -0.3387886   -0.3387886
## Mouth599    0.5184764   -0.1746708   -0.1746708

# Re-apply the CLR normalisation with a small offset
normalised.2 <- normalise.clr(Koren.16S$data.TSS, offset = 0.000001)
normalised.2[1:3, 3:5]
##                410908       177792      4294607
## Feces659   -0.3856632    4.2061195   -0.3856632
## Feces309   -0.3383694   -0.3383694   -0.3383694
## Mouth599    0.5164012   -0.1743060   -0.1743060
```

Where the intention is to apply the centered log-ratio to non-OTU data, the
function above should be avoided as it applies an inter-sample normalisation.
For this purpose OmicsFold also provides a related function
`normalise.clr.within.features()` which ensures that the means used to center
the log-ratio are calculated within a feature instead.  This can be more
appropriate for non-OTU data.
