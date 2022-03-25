#' Run BlockRank for multi-omic model
#'
#' @description
#' Takes a DIABLO model trained on multi-omic data using the mixomics package
#' and calculates the BlockRank scores for each feature.
#' This function can also calculate scores for block.plsda models.
#'
#' @param diablo.model Trained multi-omics mixomics DIABLO model
#'
#' @return A list of numeric vectors containing BlockRank scores.
#' One list item per data block
#' @export
#'
#' @examples
#' \dontrun{
#' diablo.blockrank.scores <- blockrank.diablo(diablo.model)
#' }
blockrank.diablo <- function(diablo.model) {
    y.loadings <- diablo.model[["loadings"]][["Y"]]
    x.data <- diablo.model[["X"]]
    y.data <- diablo.model[["ind.mat"]]
    ncomp <- diablo.model[["ncomp"]][1]
    nblock <- length(x.data)
    x.loadings <- diablo.model[["loadings"]][seq_len(nblock)]

    blockrank.i <-
        .internal.blockrank.calculation(x.data, y.data, x.loadings, y.loadings, ncomp, nblock)

    names(blockrank.i) <- names(x.data)

    return(blockrank.i)
}

#' Run BlockRank for single-omic model
#'
#' @description
#' Takes an sPLS-DA model trained on single-omic data using the mixomics package
#' and calculates the BlockRank scores for each feature.
#' This function can also calculate scores for plsda models.
#'
#' @param splsda.model trained single-omic mixomics sPLS-DA model
#'
#' @return A list of 1 numeric vector containing BlockRank scores. One list item for single data block
#' @export
#'
#'@examples
#' \dontrun{
#' splsda.blockrank.scores <- blockrank.splsda(splsda.model)
#' }
blockrank.splsda <- function(splsda.model) {
    x.data <- list(splsda.model[["X"]])
    y.data <- splsda.model[["ind.mat"]]
    x.loadings <- splsda.model[["loadings"]]["X"]
    y.loadings <- splsda.model[["loadings"]][["Y"]]
    ncomp <- splsda.model[["ncomp"]]
    nblock <- 1

    blockrank.i <-
        .internal.blockrank.calculation(x.data, y.data, x.loadings, y.loadings, ncomp, nblock)

    names(blockrank.i) <- c("Data")

    return(blockrank.i)
}

#' Internal BlockRank Calculation
#'
#' @param x.data List containing the centered and standardized original predictor matrix
#' @param y.data Matrix with the position of the outcome Y in the output list X
#' @param x.loadings List containing the estimated loadings for the variates.
#' @param y.loadings Matrix containing the estimated loadings for Y.
#' @param ncomp Number of components included in the model for each block
#' @param nblock Number of blocks of data
#'
#' @return A list of numeric vectors containing BlockRank scores. One list item per data block
#'
.internal.blockrank.calculation <-
    function(x.data,
             y.data,
             x.loadings,
             y.loadings,
             ncomp,
             nblock) {
        X.h <- x.data
        Y_model <- matrix(0, nrow = nrow(y.data), ncol = ncol(y.data))
        Xt_model <- lapply(x.data,
                           function(x)
                               matrix(0, nrow = ncol(x), ncol = nrow(x)))

        for (h in seq_len(ncomp)) {
            Y.b <- y.data %*% y.loadings[, h]
            Y_model <- Y_model + Y.b %*% t(y.loadings[, h])

            for (q in seq_len(nblock)) {
                X.a <- X.h[[q]] %*% x.loadings[[q]][, h]

                Xt_model[[q]] <-
                    Xt_model[[q]] + x.loadings[[q]][, h] %*% t(X.a)

                X.h[[q]] <- X.h[[q]] -
                    (X.h[[q]] %*% x.loadings[[q]][, h]) %*% t(x.loadings[[q]][, h])
            }
        }

        blockrank.i <-
            lapply(Xt_model, function(x)
                apply((x %*% Y_model), 1, function(y)
                    sum(y ^ 2)))

        blockrank.i <-
            lapply(blockrank.i, function(x)
                x / max(unlist(blockrank.i)))

        blockrank.i <- .add.feature.labels(blockrank.i, x.data)

        return(blockrank.i)
    }

#' Internal function to add feature names
#'
#' @param blockrank.i A list of numeric vectors containing BlockRank scores. One list item per data block
#' @param x.data List containing the centered and standardized original predictor matrix (with feature names)
#'
#' @return A list of numeric vectors containing BlockRank scores, with named features
#'
.add.feature.labels <- function(blockrank.i, x.data) {
    for (q in seq_along(blockrank.i)) {
        names(blockrank.i[[q]]) <- colnames(x.data[[q]])
    }
    return(blockrank.i)
}

#' BlockRank Plotting Function
#'
#' @description
#' Generates a bar plot of the top n blockrank scores, ordered highest to lowest,
#' for visualization purposes.
#' The default value 20 for nscores is not meaningful, just an example of approximately
#' how many features it is appropriate to visualise this way
#'
#' @param blockrank.i A list of numeric vectors containing BlockRank scores. One list item per data block
#' @param nscores The number of scores to display. Default 20.
#' @param feature.font.size Size of font of feature labels on y axis
#' @param model The type of model, used only in plot title.
#' @param data_source The name of the data source, used only in plot title.
#'
#' @return A ggplot2 plot object
#' @export plot.blockrank.scores
#'
#' @examples
#' \dontrun{
#' plot.blockrank.scores(blockrank.score)
#' }
plot.blockrank.scores <-
    function(blockrank.i,
             nscores = 20,
             feature.font.size = 8,
             model = "",
             data_source = "") {
        plot.data <- vector(mode = "list", length = length(blockrank.i))
        for (q in seq_along(blockrank.i)) {
            block <- names(blockrank.i)[q]
            feature <- names(blockrank.i[[q]])
            blockrank.score <- blockrank.i[[q]]
            plot.data[[q]] <-
                data.frame(block, feature, blockrank.score)
        }
        plot.data <- bind_rows(plot.data)
        plot.data <- arrange(plot.data, desc(blockrank.score))
        plot <-
            ggplot(data = head(plot.data, nscores),
                   aes(
                       x = reorder(feature, blockrank.score),
                       y = blockrank.score,
                       fill = block
                   )) +
            geom_col() +
            theme_minimal() +
            theme(
                plot.title = element_text(hjust = 0.5) ,
                plot.subtitle = element_text(hjust = 0.5),
                axis.text.y = element_text(size = feature.font.size)
            ) +
            scale_fill_brewer(palette = "Set2") +
            labs(
                x = "Feature",
                y = "BlockRank Score",
                title = sprintf("%s model of %s data", model, data_source),
                subtitle = sprintf("Top %s BlockRank Features", nscores)
            ) +
            coord_flip()

        if (length(blockrank.i) == 1) {
            plot <- plot + theme(legend.position = "none")
        }

        return(plot)
    }
