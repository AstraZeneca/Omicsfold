library("mixOmics")
library("MASS")
library(OmicsFold)
library(mixOmics)
library(dplyr)

# Singleomics case
# Create simulated model

set.seed(20220222)

n_features <- 3
n_obs <- 1e4
X <- mvrnorm(n = n_obs, mu = rep(0, n_features),
             Sigma = diag(n_features))
eps <- rnorm(n_obs, 0, 0.1)

Y <- X[, 1] + 0.5*X[, 2] + 0.3*X[, 3] + eps
Y <- as.factor(Y > 0)

splsda.simulated <- splsda(X, Y, ncomp = 2, keepX = rep(2, 2))
blockrank.splsda <- blockrank.splsda(splsda.simulated)

# Test simulated singleomics model
test_that("number of scores is number of features single", {
    expect_equal(length(unlist(blockrank.splsda)), n_features)
})

test_that("scores are between 0 and 1 single", {
    expect_true(max(unlist(blockrank.splsda)) <= 1 && min(unlist(blockrank.splsda)) >=0)
})

test_that("at least one feature has score = 1 single", {
    expect_equal(max(unlist(blockrank.splsda)), 1)
})

# Test expected ranks correct
test_that("blockrank score of x1 is greater than that of x2, which is greater than that of x3 single",{
    expect_true(blockrank.splsda[[1]][1] > blockrank.splsda[[1]][2] && blockrank.splsda[[1]][2] > blockrank.splsda[[1]][3])
})

# Create plsda model
plsda.simulated <- plsda(X, Y, ncomp = 2)
blockrank.plsda <- blockrank.splsda(plsda.simulated)

# Test simulated plsda model
test_that("number of scores is number of features plsda", {
    expect_equal(length(unlist(blockrank.plsda)), n_features)
})

test_that("scores are between 0 and 1 plsda", {
    expect_true(max(unlist(blockrank.plsda)) <= 1 && min(unlist(blockrank.plsda)) >=0)
})

test_that("at least one feature has score = 1 single", {
    expect_equal(max(unlist(blockrank.plsda)), 1)
})

test_that("blockrank score of x1 is greater than that of x2, which is greater than that of x3 plsda",{
    expect_true(blockrank.plsda[[1]][1] > blockrank.plsda[[1]][2] && blockrank.plsda[[1]][2] > blockrank.plsda[[1]][3])
})

# Multiomics tests
# Creating simulated model

X <- list(block1 = X[, 1:2], block2 = as.matrix(X[, 3]))

design <- matrix(0 , 2, 2)
ncomp <- 2

block.splsda.simulated <- block.splsda(X, Y, ncomp = ncomp, design = design)
blockrank.simulated <- blockrank.diablo(block.splsda.simulated)



# Testing simulated model

test_that("number of scores is number of features simulated", {
  expect_equal(length(unlist(blockrank.simulated)), n_features)
})

test_that("scores are between 0 and 1 simulated", {
  expect_true(max(unlist(blockrank.simulated)) <= 1 && min(unlist(blockrank.simulated)) >=0)
})

test_that("at least one feature has score = 1 simulated", {
  expect_equal(max(unlist(blockrank.simulated)), 1)
})

test_that("blockrank score of x1 is greater than that of x2, which is greater than that of x3 simulated",{
    expect_true(blockrank.simulated[[1]][1] > blockrank.simulated[[1]][2] && blockrank.simulated[[1]][2] > blockrank.simulated[[2]][1])
})

# Create block.plsda model
block.plsda.simulated <- block.splsda(X, Y, ncomp = ncomp, design = design)
blockrank.block.plsda <- blockrank.diablo(block.plsda.simulated)

# Testing block.plsda model

test_that("number of scores is number of features block.plsda", {
    expect_equal(length(unlist(blockrank.block.plsda)), n_features)
})

test_that("scores are between 0 and 1 block.plsda", {
    expect_true(max(unlist(blockrank.block.plsda)) <= 1 && min(unlist(blockrank.block.plsda)) >=0)
})

test_that("at least one feature has score = 1 block.plsda", {
    expect_equal(max(unlist(blockrank.block.plsda)), 1)
})

test_that("blockrank score of x1 is greater than that of x2, which is greater than that of x3 block.plsda",{
    expect_true(blockrank.block.plsda[[1]][1] > blockrank.block.plsda[[1]][2] && blockrank.block.plsda[[1]][2] > blockrank.block.plsda[[2]][1])
})

# Creating example model

data('breast.TCGA')
diablo.data <- list(mRNA = breast.TCGA$data.train$mrna,
                    miRNA = breast.TCGA$data.train$mirna,
                    proteomics = breast.TCGA$data.train$protein)

diablo.Y <- breast.TCGA$data.train$subtype

diablo.design <- matrix(0.1, ncol = length(diablo.data),
                        nrow = length(diablo.data),
                        dimnames = list(names(diablo.data), names(diablo.data)))
diag(diablo.design) <- 0

diablo.ncomp <- 2
diablo.keepX <- list(mRNA = c(6, 14), miRNA = c(5, 18), proteomics = c(6, 7))
example.model <- block.splsda(X = diablo.data, Y = diablo.Y,
                             ncomp = diablo.ncomp, keepX = diablo.keepX,
                             design = diablo.design)

blockrank.example <- blockrank.diablo(example.model)

# Testing example model
test_that("number of scores is number of features example", {
    n_features <- sum(sapply(example.model[["X"]], ncol))
    expect_equal(length(unlist(blockrank.example)), n_features)
})

test_that("scores are between 0 and 1 example", {
    expect_true(max(unlist(blockrank.example)) <= 1 && min(unlist(blockrank.example)) >=0)
})

test_that("at least one feature has score = 1 example", {
    expect_equal(max(unlist(blockrank.example)), 1)
})

test_that("all features with score of 0 are all the features not in model example", {
    all.features.scores <- unlist(blockrank.example)
    features.score.0 <- all.features.scores[all.features.scores == 0]
    feature.loadings <- example.model[["loadings"]][names(example.model[["loadings"]]) != "Y"]
    feature.loadings.logical <- sapply(feature.loadings, function(x) apply(x, 1, function(y) as.integer(any(as.logical(y)))))
    feature.loadings.0 <- unlist(feature.loadings.logical)
    feature.loadings.0 <- feature.loadings.0[feature.loadings.0 == 0]
    expect_equal(features.score.0, feature.loadings.0)

})


