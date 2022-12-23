#!~/.conda/envs/OmicsFold/bin/Rscript

#tune optimal number of components

args = commandArgs(trailingOnly=TRUE)
#args[1]: data
#args[2]: data labels

library(OmicsFold)
library(mixOmics)

# Read data
data<-readRDS(args[1])
data.labels<-readRDS(args[2])

#define design matrix
#values closer to 1 will identify more correlated features, but the model will be less discriminative
design <- matrix(0.1, ncol = length(data), nrow = length(data), dimnames = list(names(data), names(data)))
diag(design) = 0

# Define parameters
mfold.folds <- 10
nrepeat.find.ncomp <- 50

# Open a letter sized PDF to capture plots from this stage
pdf(file = file.path("1_find_ncomp_plots.pdf"),
    width = 11, # The width of the plots in inches
    height = 8.5) # The height of the plots in inches


# Tuning #1 ---------------------------------------------------------------
# Tune the optimum number of components up to a specified limit. Max components shouldn't be more
# than ten or so.

diablo.untrained <- mixOmics::block.splsda(X = data, Y = data.labels, ncomp = 5, design = design)
perf.untrained <- mixOmics::perf(diablo.untrained, validation = 'Mfold', folds = mfold.folds,
                       nrepeat = nrepeat.find.ncomp, progressBar = TRUE, cpus = 20)
saveRDS(perf.untrained, "1_find_ncomp_perf_result.rds")

# Pick the optimum number of components according to balanced error rate and mahalanobis distance
# If only a single component is selected, manual adjustment may be useful. Prior analysis of single omics independently may inform the optimum
# number of components to choose here.
plot(perf.untrained)

# Output the error rates for each number of components
sink("1_find_ncomp_output.txt")

cat("Dimensions of the input data:\n")
for (block.name in names(data)) {
  dimensions = dim(data[[block.name]])
  cat(sprintf("%s: %i obs x %i features\n", block.name, dimensions[1], dimensions[2]))
}

cat("\nAdvised number of components based on BER:\n")
perf.untrained$choice.ncomp$WeightedVote["Overall.BER", "mahalanobis.dist"]

sink()

# Close the PDF file
dev.off()

