.permute.data.labels <- function (lab) {
  return (sample(lab))
}

.calculate.mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

.calculate.anti.mode <- function(x) {
  uniqx <- unique(na.omit(x))
  uniqx[which.min(tabulate(match(x, uniqx)))]
}

.find.permutation.extent <- function(permuted, original) {
  num.lo <- length(original[original == levels(original)[1]])
  num.hi <- length(original[original == levels(original)[2]])

  lo <- permuted[1:num.lo]
  hi <- permuted[num.lo+1:length(permuted)]

  lo.fac <- .calculate.mode(lo)
  hi.fac <- .calculate.anti.mode(lo)

  lo.match <- length(which(lo == lo.fac))
  hi.match <- length(which(hi == hi.fac))

  return ((lo.match + hi.match) / length(permuted))
}


#' Perform a re-fit for a DIABLO model after permuting labels
#'
#' @description
#' Perform a full multi-omics model fitting and performance assessment after
#' permuting the labels associated with the classes. This can take some time.
#' This is useful when you want to assess whether a model is capable of
#' overfitting the data it was given.
#'
#' @param data Object containing input data blocks.
#' @param design DIABLO block relation design matrix.
#' @param data.labels Unpermuted data class labels.
#' @param test.keepX Array of values to test for sparse model training.
#'
#' @return Model performance balanced error rate.
#' @export
#'
train.permuted.model <- function(data, design, data.labels, test.keepX) {
  permuted.data.labels <- .permute.data.labels(data.labels)
  print ("Permuted data labels: ")
  print (permuted.data.labels)

  print ("Evaluating error rate over 6 components")
  sgccda.res.permuted <- mixOmics::block.splsda(X = data,
                                                Y = permuted.data.labels,
                                                ncomp = 6,
                                                design = design)
  perf.diablo.permuted <- mixOmics::perf(sgccda.res.permuted,
                                         validation = 'Mfold',
                                         folds = 8,
                                         nrepeat = 50,
                                         cpus = 4,
                                         progressBar = TRUE)
  ncomp = perf.diablo.permuted$choice.ncomp$WeightedVote["Overall.BER",
                                                         "centroids.dist"]
  print (paste("Optimal components: ", ncomp, sep=""))

  print ("Tuning variable penalization")
  tune.diablo.permuted <- mixOmics::tune.block.splsda(X = data,
                                                      Y = permuted.data.labels,
                                                      ncomp = ncomp,
                                                      test.keepX = test.keepX,
                                                      design = design,
                                                      validation = 'Mfold',
                                                      folds = 8,
                                                      nrepeat = 2,
                                                      cpus = 4,
                                                      dist = "centroids.dist",
                                                      progressBar = TRUE)
  list.keepX.permuted <- tune.diablo.permuted$choice.keepX
  print ("Optimal selection: ")
  print (list.keepX.permuted)

  print ("Evaluating model mfold error rate")
  sgccda.trained.permuted <- mixOmics::block.splsda(X = data,
                                                    Y = permuted.data.labels,
                                                    ncomp = ncomp,
                                                    keepX = list.keepX.permuted,
                                                    design = design)
  perf.diablo.permuted <- mixOmics::perf(sgccda.trained.permuted,
                                         validation = 'Mfold',
                                         M = 8,
                                         nrepeat = 200,
                                         dist = 'centroids.dist',
                                         cpus = 4,
                                         progressBar = TRUE)
  print ("MFold error rate: ")
  print (perf.diablo.permuted$WeightedVote.error.rate)

  return (perf.diablo.permuted$WeightedVote.error.rate)
}


#' Perform a fast model fit with permuted labels
#'
#' @description
#' Perform a quick multi-omic model fit to data with permutated class labels.
#' This does not perform the (slow) step of sparse variable or component number
#' selection - so may be more approximate.
#'
#' @param data Object containing input data blocks.
#' @param design DIABLO block relation design matrix.
#' @param data.labels Unpermuted data class labels.
#' @param ncomp Number of components to use in the model.
#' @param list.keepX.permuted Set number of variable to select from each block.
#'
#' @return List containing a representation of the permuted labels, an estimate
#' of the permutation degree (as stochastically permutation may scramble the
#' labels more or less thoroughly) and the balanced error rate over one and two
#' components.
#' @export
#'
quick.permuted.fit <- function(data, design, data.labels, ncomp,
                               list.keepX.permuted) {
  permuted.data.labels <- .permute.data.labels(data.labels)

  sgccda.trained.permuted <- mixOmics::block.splsda(X = data,
                                                    Y = permuted.data.labels,
                                                    ncomp = ncomp,
                                                    keepX = list.keepX.permuted,
                                                    design = design)
  perf.diablo.permuted <- mixOmics::perf(sgccda.trained.permuted,
                                         validation = 'Mfold',
                                         M = 8,
                                         nrepeat = 200,
                                         dist = 'centroids.dist',
                                         cpus = 4)

  error.rate.comp1 <-
      perf.diablo.permuted$WeightedVote.error.rate$centroids.dist[4,1]
  error.rate.comp2 <-
      perf.diablo.permuted$WeightedVote.error.rate$centroids.dist[4,2]

  return(list(
    "permuted.labels" = paste(as.character(permuted.data.labels),
                              collapse = " "),
    "permutation.degree" = .find.permutation.extent(permuted.data.labels),
    "error.rate.comp1" = error.rate.comp1,
    "error.rate.comp2" = error.rate.comp2
  ))
}
