#' Make a scatter plot of the degree of module membership vs the correlation with an outcome variable
#'
#' @param module Module name to plot.
#' @param module.colors List of gene module colours.
#' @param gene.module.membership Previously calculated gene module membership.
#' @param gene.trait.significance Previously calculated gene trait significance.
#' @param outcome.name Name of the outcome being correlated.
#'
#' @return Scatter plot
.plot.module.membership.correlation <- function(
    module,
    module.colors,
    gene.module.membership,
    gene.trait.significance,
    outcome.name = "outcome") {
  col.idx <- match(module, colnames(gene.module.membership))

  if (is.na(col.idx)) {
    stop("Selected module must exist in the column names of the gene.module.membership construct.")
  }

  genes.idx = which(module.colors == module)

  if (length(genes.idx) == 0) {
    stop("Selected module must exist in the list of module.colors.")
  }

  par(mfrow = c(1,1))
  WGCNA::verboseScatterplot(
      abs(gene.module.membership[genes.idx, col.idx]),
      abs(gene.trait.significance[genes.idx, 1]),
      xlab = paste("Module Membership in", module, "module"),
      ylab = paste("Gene significance for", outcome.name),
      main = paste("Module membership vs. gene significance\n"),
      cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'black'
  )
}

#' Creates a well-presented pie chart of gene module membership, with the modules represented by their actual colours
#'
#' @param moduleColors Module colours per gene.
#'
#' @return ggplot2 pie chart
#' @export
module.membership.pie <- function(moduleColors) {
  moduleCounts <- as.data.frame(table(moduleColors))

  pie.data <- moduleCounts %>%
    dplyr::mutate(end = 2 * pi * cumsum(Freq)/sum(Freq),
           start = c(0, head(end, -1)),
           middle = 0.5 * (start + end),
           hjust = ifelse(middle > pi, 1, 0),
           vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1)
    )

  pie <- ggplot2::ggplot(pie.data, ggplot2::aes(label = moduleColors))
  pie <- pie + ggforce::geom_arc_bar(ggplot2::aes(x0 = 0, y0 = 0, r0 = 0, r = 1, start = start, end = end, fill = moduleColors))
  pie <- pie + ggplot2::scale_fill_manual(values=gplots::col2hex(moduleCounts$moduleColors))
  pie <- pie + ggplot2::geom_text(ggplot2::aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = Freq, hjust = hjust, vjust = vjust)) + ggplot2::coord_fixed()
  pie <- pie + ggplot2::scale_x_continuous(limits = c(-1.4, 1.4), name = "", breaks = NULL, labels = NULL)
  pie <- pie + ggplot2::scale_y_continuous(limits = c(-1.1, 1.1), name = "", breaks = NULL, labels = NULL)
  pie <- pie + ggplot2::theme(panel.background = ggplot2::element_blank())

  return(pie)
}

#' Filter data for WGCNA analysis to remove genes with poor levels of expression or many missing values.
#'
#' @param norm.expr.data Normalised expression data for the genes, with genes in columns and samples in rows.
#' @param min.fraction The minimum fraction of samples that each gene must have data for to be considered good.
#' @param verbosity An integer level of verbosity passed to WGCNA where higher numbers create more output.
#'
#' @return A data frame containing only good samples and good genes to be analysed with the rest of WGCNA methods.
#' @export
wgcna.filter.data <- function(norm.expr.data, min.fraction = 0.75, verbosity = 1) {
  gsg <- WGCNA::goodSamplesGenes(norm.expr.data, minFraction = min.fraction, verbose = verbosity)

  cat("\nBad genes:\n")
  cat(sum(!gsg$goodGenes), "\n")

  cat("\nBad samples:\n")
  cat(sum(!gsg$goodSamples), "\n")

  return(norm.expr.data[gsg$goodSamples, gsg$goodGenes])
}

#' Perform initial clustering of genes into a dendrogram to identify outliers.
#' A horizontal red line can be added to the plot for reference where a cut point is being considered.
#'
#' @param filtered.data Filtered and normalised expression data for genes to cluster.
#' @param plot.threshold An optional height to plot a horizontal red line.
#'
#' @return The sample tree for the clustered genes.
#' @export
wgcna.initial.clustering <- function(filtered.data, plot.threshold = NULL) {
  sample.tree <- fastcluster::hclust(dist(filtered.data), method = "average")

  old.cex <- par()$cex
  par(cex = 0.6)
  plot(sample.tree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  par(cex = old.cex)

  if (!is.null(plot.threshold)) {
    abline(h = plot.threshold, col = "red")
  }

  return(sample.tree)
}

#' Cut a sample tree at a specified height to remove outliers from an initial clustering.
#'
#' @param sample.tree A sample tree to cut outliers from.
#' @param cut.threshold The height in the tree to perform the cut.
#' @param filtered.data The original expression data represented by the sample tree.
#'
#' @return Expression data containing the largest subset of genes from the sample tree after performing the cut.
#' @export
wgcna.cut.tree <- function(sample.tree, cut.threshold, filtered.data) {
  clust <- WGCNA::cutreeStatic(sample.tree, cutHeight = cut.threshold, minSize = 10)
  largest.cluster <- names(sort(table(clust), decreasing = TRUE))[1]

  return(filtered.data[(clust == largest.cluster),,drop=F])
}

#' Find the soft threshold for forming the network of modules.
#'
#' @param expr.data Filtered and normalised expression data for genes to cluster.
#' @param powers A vector of powers to trial when picking a soft threshold.
#' @param block.size The size of blocks to perform the analysis with.
#'                   NULL indicates that the system should pick a block size.
#'                   It is better to process all genes simultaneously if memory allows.
#' @param marked.threshold An optional value to mark a horizontal red line on the output plot.
#'
#' @details A plot is generated showing the scale free topology fit index for each provided power value.
#'          The highest peak at the lowest power value should be selected for network construction.
#'          Ideally this will be above 0.90.
#'
#' @export
wgcna.fit.soft.threshold <- function(
    expr.data,
    powers = c(1:10, seq(12, 20, 2)),
    block.size = NULL,
    marked.threshold = 0.90) {
  sft <- WGCNA::pickSoftThreshold(expr.data, powerVector = powers, verbose = 5, corFnc = WGCNA::bicor, blockSize = block.size)

  par(mfrow = c(1,2))
  cex1 = 0.9

  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")

  if (is.numeric(marked.threshold)) {
    abline(h = marked.threshold, col="red")
  }

  plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

  par(mfrow = c(1,1))
}

#' Construct a WGCNA network on the expression data using a specified power and block size.
#'
#' @param expr.data Filtered and normalised expression data for genes to cluster.
#' @param power The soft-thresholding power to use. This is usually fitted using `wgcna.fit.soft.threshold`.
#' @param max.block.size The maximum size for a block during network construction.
#'                       Ideally this will equal the number of genes, but will be restricted by RAM availability.
#' @param TOM.file.path An optional file path to save the TOM file(s) to. One file will be generated per block.
#' @param verbosity An integer indicating the verbosity of the WGCNA library when creating the network.
#'
#' @return The constructed network object.
#' @export
wgcna.create.blockwise.network <- function(expr.data, power, max.block.size = 6000, TOM.file.path = NULL, verbosity = 1) {
  save.TOMs <- FALSE
  save.TOM.file.base <- 'blockwiseTOM'

  if (!is.null(TOM.file.path)) {
    save.TOMs <- TRUE
    save.TOM.file.base <- file.path(TOM.file.path, 'blockwiseTOM')
  }

  return(WGCNA::blockwiseModules(
    expr.data,
    power = power,
    maxBlockSize = max.block.size,
    TOMType = "unsigned",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    saveTOMs = save.TOMs,
    saveTOMFileBase = save.TOM.file.base,
    corType = "bicor",
    pearsonFallback = "individual",
    verbose = verbosity
  ))
}

#' Extract module colors from a constructed network.
#'
#' @param network A constructed network to extract module colors from.
#'
#' @return A vector of color names equal in length and order to the genes in the network.
#' @export
wgcna.module.colors <- function(network) {
  WGCNA::labels2colors(network$colors)
}

#' Plot a dendrogram for a constructed network with module colors plotted as a horizontal stripe beneath.
#'
#' @param network A constructed network to plot a dendrogram from.
#'
#' @export
wgcna.plot.network.dendro <- function(network) {
  mergedColors <- wgcna.module.colors(network)
  WGCNA::plotDendroAndColors(
    network$dendrograms[[1]], mergedColors[network$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05
  )
}

#' Plot a pie chart of module sizes given the module colors from a constructed network.
#'
#' @param module.colors A vector or colors matching the genes from a constructed network.
#'
#' @return A plot object to render the pie chart with.
#' @export
wgcna.module.pie.plot <- function(module.colors) {
  return(module.membership.pie(module.colors))
}

#' Identify eigengenes that represent each module of genes.
#'
#' @param expr.data Filtered and normalised expression data for genes to calculate eigengenes for.
#' @param module.colors A vector of color names identifying modules of genes, the length of which matches the number of columns in `expr.data`.
#'
#' @return The eigengenes for each module, one per input sample.
#' @export
wgcna.determine.module.eigengenes <- function(expr.data, module.colors) {
  MEs0 <- WGCNA::moduleEigengenes(expr.data, module.colors)$eigengenes
  MEs <- WGCNA::orderMEs(MEs0)
  colnames(MEs) <- substring(colnames(MEs), 3)

  return(MEs)
}

#' Calculate the biweight midcorrelation between two matrices of data for the same set of samples.
#'
#' @param x Matrix x for the correlation.
#' @param y Matrix y for the correlation.
#'
#' @return A list containing the correlation values and the p values for each correlation between the two matrices.
#' @export
wgcna.calculate.correlation <- function(x, y) {
  correlation <- WGCNA::bicor(x, y, use = "p")
  p.values <- WGCNA::bicorAndPvalue(x, y, use = "p")$p

  return(list(correlation = correlation, p.values = p.values))
}

.var.names.to.title.case <- function(var.names) {
  var.names %>%
    stringr::str_replace('\\.', ' ') %>%
    tools::toTitleCase()
}

#' Plot a labelled heatmap of correlations, sorted by correlation value.
#'
#' @param correlation A list containing the correlation values and p values for the correlations.
#'
#' @export
wgcna.plot.labelled.heatmap <- function(correlation) {
  data.matrix <- as.matrix(correlation$correlation)
  y.labels <- rownames(correlation$correlation)

  # Make x.labels be more presentable
  x.labels <- colnames(correlation$correlation) %>%
    .var.names.to.title.case()

  # Workaround bug where WGCNA doesn't handle single column input matrix for heatmap
  if (ncol(data.matrix) == 1) {
    data.matrix <- cbind(data.matrix, data.matrix)
    x.labels <- c(x.labels, x.labels)
  }

  text.matrix <- sprintf("%0.2f  (p = %g)", signif(data.matrix, 2), signif(correlation$p.values, 1))
  dim(text.matrix) <- dim(data.matrix)

  # Sort the modules by correlation value for a cleaner plot
  sort.order <- order(correlation$correlation[,1], decreasing=TRUE)
  y.labels <- y.labels[sort.order]
  data.matrix <- data.matrix[sort.order,]
  text.matrix <- text.matrix[sort.order,]

  WGCNA::labeledHeatmap(
    Matrix = data.matrix,
    xLabels = x.labels,
    yLabels = y.labels,
    ySymbols = y.labels,
    colorLabels = FALSE,
    colors = rev(WGCNA::greenWhiteRed(50)),
    textMatrix = text.matrix,
    setStdMargins = FALSE,
    cex.text = 1.0,
    zlim = c(-1,1),
    main = "Module-trait relationships"
  )
}

#' Generate all the plots of module membership correlation for all modules.
#'
#' @param module.colors A list of module color names the same length as the number of gene expression data in the correlations.
#' @param expr.module.correlation A correlation between the gene expression data and the module eigengenes.
#' @param expr.trait.correlation A correlation between the gene expression data and a trait from the sample metadata.
#' @param include.modules A list of module to include in the plots, or NULL for all modules.
#' @param ignore.modules A list of module names to ignore when creating plots.
#'
#' @return A GO enrichment object which contains objects such as the `enrichmentTable` which can be exported for later analysis of associated GO terms for each module.
#' @export
wgcna.plot.module.membership.correlations <- function(
    module.colors,
    expr.module.correlation,
    expr.trait.correlation,
    include.modules = NULL,
    ignore.modules = c('grey')) {
  if (is.null(include.modules)) {
    include.modules <- colnames(expr.module.correlation$correlation)
  }

  for (module.name in include.modules) {
    if (module.name %in% ignore.modules) {
      next
    }

    .plot.module.membership.correlation(
      module.name,
      module.colors,
      expr.module.correlation$correlation,
      expr.trait.correlation$correlation,
      outcome = .var.names.to.title.case(colnames(expr.trait.correlation$correlation)))
  }
}

#' Provide GO term enrichment for modules identified in the network construction.
#'
#' @param gene.list A vector listing the names of genes included in modules to find GO terms for.
#' @param module.labels A vector of labels the same length as `gene.list` identifying the module each gene belongs to.
#' @param organism The name of the source organism (e.g. mouse).
#' @param sub.collections An optional list of sub-collections for GO terms to be included in the analysis (e.g. GO.BP for biological pathways).
#' @param threshold.type The thershold type to use for multiple testing p-value correction. FDR is the recommended type, but bonferroni can also be chosen.
#' @param ignore.labels A list of module labels to ignore during the analysis.
#' @param best.data.sets.count The number of best GO terms to report in the output for each module.
#'
#' @return A GO enrichment object which contains objects such as the `enrichmentTable` which can be exported for later analysis of associated GO terms for each module.
#' @export
enrich.modules <- function(
    gene.list,
    module.labels,
    organism = "mouse",
    sub.collections = c("GO.BP"),
    threshold.type = "FDR",
    ignore.labels = c("grey"),
    best.data.sets.count = 10) {
  # Convert gene names to entrez IDs
  gene.symbols <- map.symbols.entrez.ids(gene.list, organism = organism)

  # Prepare our reference collection of GO terms
  GO.collection <- anRichment::buildGOcollection(organism = organism)
  if (!is.null(sub.collections) & length(sub.collections) > 1) {
    GO.collection <- anRichmentMethods::subsetCollection(GO.collection, tags = sub.collections)
  }

  # Enrich the list of genes
  GO.enrichment = anRichmentMethods::enrichmentAnalysis(
    classLabels = module.labels,
    identifiers = gene.symbols$entrez.id,
    refCollection = GO.collection,
    useBackground = "given",
    threshold = 1e-4,
    thresholdType = threshold.type,
    getOverlapEntrez = TRUE,
    getOverlapSymbols = TRUE,
    ignoreLabels = ignore.labels,
    nBestDataSets = best.data.sets.count
  )

  return(GO.enrichment)
}
