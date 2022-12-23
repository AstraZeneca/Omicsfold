#!~/.conda/envs/OmicsFold/bin/Rscript

#tune optimal number of features per component

args = commandArgs(trailingOnly=TRUE)
#args[1]: data
#args[2]: data labels
#args[3]: perf untrained

library(OmicsFold)
library(mixOmics)
library(BiocParallel)
library(batchtools)

# Read data
data<-readRDS(args[1])
data.labels<-readRDS(args[2])
perf.untrained <- readRDS(args[3])

#define design matrix
#values closer to 1 will identify more correlated features, but the model will be less discriminative
design <- matrix(0.1, ncol = length(data), nrow = length(data), dimnames = list(names(data), names(data)))
diag(design) = 0

# Define parameters

mfold.folds <- 10
nrepeat.keepX <- 50

#determine optimal number of components to mahalanobis distance
ncomp.keepX <- perf.untrained$choice.ncomp$WeightedVote["Overall.BER", "mahalanobis.dist"]

#overwrite to ncomp to 2 if 1
if (ncomp.keepX == 1){
  ncomp.keepX = 2
}

# Set to test a range of variables to test
#test 10%, 20%, 30%, 40% of features
keepx_grid <- function(i) {
 block_length <- ncol(data[[names(data[i])]])
 list.keepX <- c(round(block_length/10, digits=0),  round(block_length/5, digits=0),  round(block_length/3.3, digits=0), round(block_length/2.5, digits=0))
 return(list.keepX)
}

test.keepX <- lapply(names(data),keepx_grid)
names(test.keepX) <- names(data)

#loop over each block, and if number of features is greater than max available, use 50% of available features
for (i in names(test.keepX)) {
  max_features <- ncol(data[[i]])
  for(x in 1:length(test.keepX[[i]])){
    print(test.keepX[[i]][[x]])
    if (test.keepX[[i]][[x]] > max_features){
    test.keepX[[i]][[x]] <- round(max_features/2)
    }
  }  
}


# Open a letter sized PDF to capture plots from this stage
pdf(file = "2_find_keepX_plots.pdf",
    width = 11, # The width of the plots in inches
    height = 8.5) # The height of the plots in inches


# Tuning #2 ---------------------------------------------------------------

# Establish the optimum number of parameters to retain in the sparse model. DIABLO uses lasso penalization to select variables, but the mixOmics
# implementation allows the user to specify a number of variables rather than lamdba parameters (which can be difficult to estimate).
# The optimum number of parameters to retain has similar considerations to sPLS-DA - aim for 1/3 of total variables, and ensure the low range is well covered
# to allow elimination of noisy variable as much as possible.
# Every component of every block can have a different penalization parameter, so this step can be very computationally intensive to run. As for component
# numbers, prior analysis of single omics can inform the range of keepX values to test here and so cut down on the number that have to be processed.

if (!exists("already.tested.X")) {
  already.tested.X <- list()
  for (block.name in names(data)) {
    already.tested.X[[block.name]] <- vector()
  }
}


tune.diablo <- tune.block.splsda(X = data, Y = data.labels, ncomp = ncomp.keepX, test.keepX = test.keepX, design = design,
                                 validation = 'Mfold', folds = mfold.folds, nrepeat = nrepeat.keepX,
                                dist = "mahalanobis.dist", progressBar = FALSE, BPPARAM = BatchtoolsParam(workers = 20, cluster="slurm"))

saveRDS(tune.diablo, "2_find_keepX_tune_result.rds")

# Find the error rate for the tuning at different parameter values, the optimal number of components, and the number of selected variables for each component
sink("2_find_keepX_output.txt")

cat("\n\n")
cat("keepX values tested:\n")
test.keepX

cat("\n\n")
cat("Advised counts for keepX:\n")
tune.diablo$choice.keepX

ber.per.component <- function(keepX.results) {
  ber.per.component <- vector()

  last.comp <- length(keepX.results$choice.keepX[[1]])
  first.comp <- last.comp - ncol(keepX.results$error.rate) + 1

  for (comp in first.comp:last.comp) {
    keepX.values <- array()
    for (block.index in 1:length(keepX.results$choice.keepX)) {
      keepX.values[block.index] <- keepX.results$choice.keepX[[block.index]][comp]
    }

    row.key <- paste(keepX.values, collapse="_")
    col.key <- sprintf("comp%i", comp)

    ber.per.component[col.key] <- keepX.results$error.rate[row.key, col.key]
  }

  return(ber.per.component)
}

cat("\n")
cat("Error rates for selected keepX values:\n")
print(ber.per.component(tune.diablo))

cat("\n")
cat("Advised number of components:\n")
print(tune.diablo$choice.ncomp)

cat("Error rates per component / keepX value combination:\n")
tune.diablo$error.rate

sink()

plot(tune.diablo, col = color.jet(ncomp.keepX))

# Close plots PDF
dev.off()

