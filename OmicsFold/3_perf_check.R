#!~/.conda/envs/OmicsFold/bin/Rscript

#Determine model performance using Mfold cross-validation 
#generate model visualizations 

args = commandArgs(trailingOnly=TRUE)
#args[1]: data
#args[2]: data labels
#args[3]: tune diablo

library(OmicsFold)
library(mixOmics)
library(dplyr)

# Read data
data<-readRDS(args[1])
data.labels<-readRDS(args[2])
tune.diablo<-readRDS(args[3])

#define design matrix
#values closer to 1 will identify more correlated features, but the model will be less discriminative
design <- matrix(0.1, ncol = length(data), nrow = length(data), dimnames = list(names(data), names(data)))
diag(design) = 0

# Define parameters
mfold.folds <- 10
nrepeat.tuned <- 50

ncomp.tuned <- tune.diablo$choice.ncomp$ncomp

#overwrite to ncomp to 2 if 1

if (ncomp.tuned == 1){
  ncomp.tuned = 2
}

keepX.tuned <- tune.diablo$choice.keepX

#make sure keepx.tuned is not greater in length than ncomp.tuned
if (length(keepX.tuned[[1]]) != ncomp.tuned){
  for (i in 1:length(keepX.tuned)){
    keepX.tuned[[i]] <-  keepX.tuned[[i]][1:ncomp.tuned]
  }
}

save.core.diablo.stats <- function(file, trained.model, perf.result, centroids) {
  sink(file)
  
  cat("Final keepX:\n")
  print(keepX.tuned)

  cat("Overall error rates:\n")
  print(perf.result$WeightedVote.error.rate$mahalanobis.dist)

  cat("\n")
  cat("Error rates per block and component during DIABLO performance test:\n")
  print(perf.result$error.rate)

  for (block.name in names(trained.model$X)) {
    if (trained.model$ncomp[block.name] > 1) {
      cat(sprintf("Model variance for block '%s':\n", block.name))
      print(get.model.variance(trained.model, block = block.name))
      cat("\n")
    }
  }

  cat("Centroids consensus table:\n")
  print(centroids)

  sink()
}

save.loadings.tables <- function(file, trained.model, perf.result) {
  loadings.final <- data.frame()
  for (block.name in names(trained.model$X)) {
    loadings.all <- get.diablo.top.loadings.with.stability(trained.model, perf.result, block.name, feature.count = Inf)
    loadings.all$block <- block.name
    loadings.final <- rbind(loadings.final, loadings.all)
  }
    write.csv(loadings.final, file = "3_loadings_stability_all.csv")
  }


save.diablo.model.plots <- function(file, trained.model, perf.result, centroids) {
  # Open a letter sized PDF to capture plots
  pdf(file = file,
      width = 11, # The width of the plots in inches
      height = 8.5) # The height of the plots in inches

  ncomp.tuned <- trained.model$ncomp[1]

  # Component comparisons
  if (ncomp.tuned > 1) {
    for (comp1 in 1:(ncomp.tuned - 1)) {
      for (comp2 in (comp1 + 1):ncomp.tuned) {
        comps <- c(comp1, comp2)

        # Individual block sample plot
        plotIndiv (trained.model, comp = comps, group = trained.model$Y, ind.names = FALSE, ellipse = TRUE,
                  legend = TRUE, title = 'DIABLO')

        # The arrow plot, in this case, has multiple tips corresponding to each block
        plotArrow(trained.model, comp = comps, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
      }
    }
  }
  # Plot the consensus plot
  ggplot2::ggplot(centroids,
                  ggplot2::aes(comp1, comp2, color = label))+ 
    ggplot2::geom_point() + ggplot2::stat_ellipse()

  # Block correlation plot for each component
  for (comp in 1:ncomp.tuned) {
    # Unfortunately it's not possible to include a title on these plots
    plotDiablo(trained.model, ncomp = comp)
  }

  # Loadings. Specify the component to examine loadings from variables selected across each block
  for (block.name in names(trained.model$X)) {
    for (comp in 1:ncomp.tuned) {
      name.var <- Map(center.truncate, colnames(trained.model$X[[block.name]]))
      plotLoadings(trained.model, block = block.name, comp = comp, contrib = 'max', method = 'median',
                   col.ties = "grey", name.var = name.var, ndisplay=20,
                   title = sprintf("%s, comp %i", block.name, comp))
    }
  }

  # Plot feature stability for a summary overview
  for (block.name in names(trained.model$X)) {
    for (comp in 1:ncomp.tuned) {
      plot(diablo.selection.stability(perf.result, comp = comp, block=block.name), type = 'h',
          ylab = 'Stability', xlab = 'Features', las = 2,
          main = sprintf("Selection stability, block = '%s', Comp %i", block.name, comp))
    }
  }

  # Features against components
  if (ncomp.tuned > 1) {
    for (comp1 in 1:(ncomp.tuned - 1)) {
      for (comp2 in (comp1 + 1):ncomp.tuned) {
        plotVar(trained.model, var.names = FALSE, style = 'graphics', legend = FALSE, overlap = FALSE, comp = c(comp1, comp2))
      }
    }
  }

  # Circular correlation plot between blocks
  circosPlot(trained.model, cutoff = 0.6, line = TRUE, size.labels = 1, comp=1, size.variables=0.0001)

  # Clustered image map, showing multi-omics clustering
  cimDiablo(trained.model, size.legend = 0.7, comp=1, color.Y = rainbow(nlevels(trained.model$Y)))
  
  dev.off()
}


# Create model objects ----------------------------------------------------
# Train the model. This step is quick!

sgccda.trained <- block.splsda(X = data, Y = data.labels, ncomp = ncomp.tuned,
                                                  keepX = keepX.tuned, design = design)

saveRDS(sgccda.trained, "3_trained_model_object.rds")

# Determine model performance. This uses Mfold cross-validation - M will generally have to be kept quite small for the low number of samples. nrepeat is
# bound only by computational time

perf.diablo <- perf(sgccda.trained, validation = 'Mfold', folds = mfold.folds,
                   dist = 'mahalanobis.dist', nrepeat = nrepeat.tuned, progressBar = TRUE, cpus = 20)

saveRDS(perf.diablo, "3_perf_result_object.rds")

# Centroids from all blocks
#modified from get.block.centroids because function does not work with updated version of mixomics
arrow <-  mixOmics::plotArrow(sgccda.trained)
centroids <- arrow[["data"]]
centroids$sample <- rownames(centroids)
centroids <- centroids %>% select(sample, group, x_centroid, y_centroid)
colnames(centroids) <- c("sample", "label", "comp1", "comp2")
rownames(centroids) <- NULL

save.core.diablo.stats("3_perf_check_output.txt",
                       trained.model = sgccda.trained,
                       perf.result = perf.diablo,
                       centroids = centroids)

save.loadings.tables(file = "3_loadings_stability_all.csv",
                     trained.model = sgccda.trained,
                     perf.result = perf.diablo)

save.diablo.model.plots(file ="3_perf_check_plots.pdf",
                        trained.model = sgccda.trained,
                        perf.result = perf.diablo,
                        centroids = centroids)

pdf(file = "3_blockrank.pdf",
    width = 11, # The width of the plots in inches
    height = 8.5) # The height of the plots in inches

#blockrank
plot.blockrank.scores(blockrank.diablo(sgccda.trained), feature.font.size = 12)

dev.off()

