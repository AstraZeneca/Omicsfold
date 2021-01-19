
#' Internal function for dividing over complete samples in Total Sum Scaling
#' (TSS).
#'
#' @param x Row to divide over.
#'
#' @return TSS scaled row.
.TSS.divide = function(x){
	if (sum(x) > 0)
		return (x/sum(x))
	else
		return (x)
}


#' Prefilter omics analysis input data in count form (e.g. OTUs) to remove
#' features which have a total count less than a (small) proportion of the total
#' measured counts. The default threshold is one part in 10,000 (0.01\%) - this
#' is usually sufficient to remove very low-count variables, which will be
#' unreliable features for model prediction.
#'
#' @param otu.counts OTU count data frame of size n (sample) x p (OTU).
#' @param percent Cutoff chosen in percent, default to 0.01.
#'
#' @return Data frame of input data, filtered to omit features below the count
#' proportion threshold.
#' @export
#'
#' @examples
#' \dontrun{
#' low.count.filter(raw.count)
#' }
low.count.removal = function(otu.counts, percent=0.01 ) {
	keep.otu <- which(colSums(otu.counts)*100/(sum(colSums(otu.counts))) > percent)
	data.filter <- otu.counts[,keep.otu]
	return(data.filter)
}


#' Apply Total Sum Scaling (TSS) normalisation to count data, to account for
#' differences in count (e.g. sequencing) depths between samples. Giving
#' proportion of total sample counts, this is the conventional way of
#' normalising OTU count data. Optionally include an offset to avoid
#' division/log zero problems - this defaults to zero, but 1 (count) is usually
#' appropriate for any count data with totals of thousands of counts or more.
#'
#' @param otu.counts OTU count data frame of size n (sample) x p (OTU).
#' @param offset Offset to apply, defaulting to zero.
#'
#' @return Data frame containing count data normalised as proportion of total
#' sample counts.
#' @export
#'
#' @examples
#' \dontrun{
#' normalise.tss(otu.count, offset=1)
#' }
normalise.tss = function(otu.counts, offset=0) {
	offset <- otu.counts + offset
	return(t(apply(offset, 1, .TSS.divide)))
}


#' Alternate Cumulative Sum Scale (CSS) method for normalising count data for
#' inter-sample depth. Relies upon the metagenomeSeq implementation.
#'
#' @param otu.counts OTU count data frame of size n (sample) x p (OTU).
#'
#' @return Data frame containing count data normalised cumulatively.
#' @export
#'
#' @examples
#' \dontrun{
#' normalise.css(otu.count)
#' }
normalise.css = function(otu.counts) {
	data.metagenomeSeq <- metagenomeSeq::newMRexperiment(t(otu.counts), featureData=NULL, libSize=NULL, normFactors=NULL)
	p <- metagenomeSeq::cumNormStat(data.metagenomeSeq)
	data.cumnorm <- metagenomeSeq::cumNorm(data.metagenomeSeq, p=p)
	otu.css <- t(metagenomeSeq::MRcounts(data.cumnorm, norm=TRUE, log=TRUE))
	return(otu.css)
}

#' Internal function for "empirical" logit normalisation of a feature (column)
#' of data. The empirical logit function differs for standard logit
#' normalisation in that an epsilon factor is added to ensure that function does
#' not tend to +/- infinity for input values close to 100\% and 0\%
#' respectively.
#'
#' @param feature Feature column.
#'
#' @return Normalised feature column.
.normalise.logit.feature = function(feature) {
	epsilon.min <- min(feature)
	epsilon.max <- 1-max(feature)
	epsilon <- min(epsilon.min, epsilon.max)

	# Set minimum and maximum values for the smoothing factor epsilon
	epsilon <- max(epsilon, 0.01)
	epsilon <- min(epsilon, 0.1)

	return(log((feature + epsilon)/(1 - feature + epsilon)))
}

#' Apply the empirical logit normalisation to a data frame of omics input data.
#' This is intended to convert compositional data, e.g. proportional data in the
#' range 0..1, to Euclidean space which is most appropriate for the linear
#' models. The empirical logit function differs for standard logit normalisation
#' in that an epsilon factor is added to ensure that function does not tend to
#' +/- infinity for input values close to 100\% and 0\% respectively. The logit
#' or empirical logit function will be a more appropriate choice than centred
#' log-ratio (CLR) for non-OTU data.
#'
#' @param input Data frame of input compositional data to normalise. Input data
#' should be proportions 0-1.
#'
#' @return Data normalised using empirical logit. Proportions below 0.5 will be
#' negative, but output will not tend to infinity for zero or 1 input.
#' @export
#'
#' @examples
#' \dontrun{
#' normalise.logit.empirical(data.proportional)
#' }
normalise.logit.empirical = function(input) {
	normalised <- apply(input, 2, .normalise.logit.feature)
	rownames(normalised) <- rownames(input)
	return(normalised)
}

#' Apply the standard logit normalisation to a data frame of omics input data.
#' This is intended to convert compositional data, e.g. proportional data in the
#' range 0..1, to Euclidean space which is most appropriate for the linear
#' models. The logit function will tend to +/- infinity for input values close
#' to 100% and 0% respectively. The logit or empirical logit function will be a
#' more appropriate choice than centred log-ratio (CLR) for non-OTU data.
#'
#' @param input Data frame of input compositional data to normalise. Input data
#' should be proportional 0-1.
#'
#' @return Data normalised using empirical logit. Proportions below 0.5 will be
#' negative, and output will tend to -/+ infinity for zero or 1 input.
#' @export
#'
#' @examples
#' \dontrun{
#' normalise.logit(data.proportional)
#' }
normalise.logit = function(input) {
	return(log(input/(1 - input)))
}


#' Apply centred log-ratio (CLR) normalisation to sum scaled OTU count data.
#' This is another method for converting the compositional data, i.e.
#' proportional data in the range 0..1 to Euclidean space which is most
#' appropriate for the linear models. Note that this should only be applied to
#' OTU data, as it applies another inter-sample normalisation.
#'
#' @param input Scaled OTU data as proportions 0-1, e.g. output by
#' normalise.TSS().
#' @param offset Optional offset to apply to raw data to avoid logging of zero
#' values. Only needed if any zeroes are present - should generally be set very
#' small, e.g. 0.000001.
#'
#' @return Data normalised by the CLR method.
#' @export
#'
#' @examples
#' \dontrun{
#' normalise.clr(otu.data.tss)
#' }
normalise.clr = function(input, offset=0) {
  normalised.clr <- mixOmics::logratio.transfo(X = as.matrix(input), logratio = 'CLR', offset=offset)
  # Annoyingly, the output object does not allow direct access the matrix of results. This is an easy way to return it.
  return(normalised.clr[,])
}

#' Apply centred log-ratio (CLR) normalisation to other compositional data, but
#' restrict normalisation to *within* features only. This is another method for
#' converting the compositional data, i.e. proportional data in the range 0..1
#' to Euclidean space which is most appropriate for the linear models. The
#' implementation is the same as CLR, but on transposed input data (which is
#' then transposed back to the input orientation). Note this is experimental,
#' though it does give a sensible normalisation.
#'
#' @param input Data as proportions 0-1.
#' @param offset Optional offset to apply to raw data to avoid logging of zero
#' values. Only needed if any zeroes are present - should generally be set very
#' small, e.g. 0.000001.
#'
#' @return Data normalised by the within-feature CLR method
#' @export
#'
#' @examples
#' \dontrun{
#' normalise.clr.within.features(data.proportional)
#' }
normalise.clr.within.features = function(input, offset=0) {
  normalised.clr <- mixOmics::logratio.transfo(X = t(as.matrix(input)), logratio = 'CLR', offset=offset)
  return(t(normalised.clr[,]))
}
