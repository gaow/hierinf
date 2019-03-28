require("testthat")

## random number generator
RNGkind("L'Ecuyer-CMRG")

test_that("test_hierarchy: check input", {
  expect_error(test_hierarchy(x = NULL, y = NULL, dendr = NULL, family = NULL),
               "The response y is required to be a vector or a list of vectors if multiple data sets are present.")

  expect_error(test_hierarchy(x = NULL, y = matrix(1:4, ncol = 2), dendr = NULL,
                              family = NULL),
               "The elements of the list y are required to be numeric vectors or matrices with only one column. In the case of only one data set, it is enough that y is a numeric vector or matrix with only one column but it can as well be a list with one element.")

  expect_error(test_hierarchy(x = NULL, y = list(matrix(1:4, ncol = 2),
                                                 matrix(1:4, ncol = 2)),
                              dendr = NULL, family = NULL),
               "The elements of the list y are required to be numeric vectors or matrices with only one column.")

  expect_error(test_hierarchy(x = NULL, y = 1:2, dendr = NULL,
                              family = NULL),
               "The input x is required to be a matrix or a list of matrices if multiple data sets are present.")

  expect_error(test_hierarchy(x = matrix(1:4, ncol = 2), y = 1:2, dendr = NULL,
                              family = NULL),
               "The matrix x is required to have column names. If there is no natural naming convention, then one can set them to some integer, say, 1 to p.",
               fixed = TRUE)

  tt <- matrix(1:4, ncol = 2)
  colnames(tt) <- c("a", "b")
  expect_error(test_hierarchy(x = tt, y = 1:2, dendr = NULL,
                              family = NULL),
               "The input dendr is required to be a list of dendrograms.")

  # Please note that family = NULL results in taking the default value.

  set.seed(88)
  x <- matrix(rnorm(100), ncol = 2)
  colnames(x) <- c("col1", "col2")
  y <- 1:50
  dendr <- cluster_var(x = x)
  # expect_error(test_hierarchy(x = x, y = y, dendr = dendr, family = "alsdkk"),
  #              "'arg' should be one of \"gaussian\", \"binomial\"")

  set.seed(124)
  res.multisplit <- multisplit(x = x, y = y)
  expected_result_1 <- data.frame(block = NA, p.value = NA,
                                  significant.cluster = NA)
  expected_result_1$significant.cluster <- list(NA)
  attr(expected_result_1, "class") <- c("data.frame")
  expected_result <- list(res.multisplit = res.multisplit,
                          res.hierarchy = expected_result_1)
  attr(expected_result, "class") <- c("hierT", "list")
  expect_equal({set.seed(124); test_hierarchy(x = x, y = y, dendr = dendr,
                                              family = "gaussian")},
               expected_result)
})

#### Check output with one data set ####

# This function calculates the p-value for each of the notes in the tree.
check_test_hierarchy <- function(x, y, clvar, res.multisplit, B, cluster_test){
  CT_colnames <- lapply(X = cluster_test, function(x) paste(x, collapse = "_"))
  res <- matrix(NA, ncol = length(cluster_test), nrow = B)
  colnames(res) <- CT_colnames
  RES <- rep(NA, length(cluster_test))
  names(RES) <- CT_colnames
  for (i in cluster_test) { # for each cluster
    for (b in 1:B) { # for each split (multi-sample splitting)
      # selected coefficients
      sel.coef <- res.multisplit[[1]]$sel.coef[b, ][!is.na(res.multisplit[[1]]$sel.coef[b, ])]
      # other half of the samples
      ind <- res.multisplit[[1]]$out.sample[b, ]
      # combined data set (clvar plus x)
      clvar_x <- cbind(clvar, x[, sel.coef, drop = FALSE])[ind, ]

      # intersection and set difference of selected coefficients and the given cluster
      intersect_i <- intersect(sel.coef, i)
      setdiff_i <- setdiff(sel.coef, i)

      # columnnames for the model \hat{S}^{(b)} \setminus C
      sel_i_clvar <- c(colnames(clvar), setdiff_i)

      # browser()
      # design matrix of the reduced model: variables of \hat{S}^{(b)} \setminus C & clvar
      clvar_x_reduced <- clvar_x[, sel_i_clvar, drop =  FALSE]
      if (ncol(clvar_x_reduced) == 0) {
        clvar_x_reduced <- rep(1, length(y[ind]))
      }

      res[b, paste(i, collapse = "_")] <-
        if (length(intersect_i) == 0) { # Equation (2) on page 333 of Mandozzi and Buehlmann (2016)
          1
        } else {
          # min(1, ... * |\hat{S}^{(b)}| / |\hat{S}^{(b)} \setminus C|)  =>  Equation (2) & (3)
          # on page 333 of Mandozzi and Buehlmann (2016)
          min(1, anova(
            # full model: variables of \hat{S}^{(b)} & clvar
            lm(y[ind] ~ clvar_x),
            # reduced model: variables of \hat{S}^{(b)} \setminus C & clvar
            lm(y[ind] ~ clvar_x_reduced),
            test = "F")$P[2] *
              length(sel.coef) / length(intersect_i)
          )
        }
    }
    # Equation (4) on page 333 of Mandozzi and Buehlmann (2016)
    RES[paste(i, collapse = "_")] <- adj_pval(res[,  paste(i, collapse = "_")], B = B)
  }

  # hierarchical adjustment has to be done by hand (Equation below Equation (4))
  return(RES)
}

adj_pval <- function(pvals, B) {
  # define the sequence of gamma values
  gamma_min <- 0.05
  gamma_step <- 0.01
  gamma_seq <- seq(gamma_min, 1, gamma_step)

  # compute the empirical quantile vector
  gamma_step <- vector("numeric", length = length(gamma_seq))
  for (g in 1:length(gamma_seq)) {
    gamma_step[g] <- min(1, quantile(pvals / gamma_seq[g], gamma_seq[g],
                                     na.rm = TRUE))
  }

  # compute the adjusted p value
  # Equation 4 on page 333 in Mandozzi and Buehlmann (2016)
  return(min(1, (1 - log(gamma_min)) * min(gamma_step)))
}

### Example I ###
test_that("test_hierarchy: check output (Example I)", {
  ## simulate index
  n <- 800
  p <- 5
  B <- 50

  ## simulate data
  require(MASS)
  set.seed(9229)
  sim.geno <- mvrnorm(n = n, mu = rep(0, p),
                      Sigma = toeplitz(0.8^(seq(0, p - 1))))
  colnames(sim.geno) <- paste0("rsid", 1:p)

  set.seed(144)
  data.dim <- dim(sim.geno)

  ind.active <- sample(1:data.dim[2], 2)
  beta <- rep(0, data.dim[2])
  beta[ind.active] <- 2
  y <- sim.geno %*% beta + rnorm(data.dim[1])

  # cluster the data
  dendr <- cluster_var(x = sim.geno)
  # plot(dendr$res.tree[[1]])

  # multisplit
  set.seed(2)
  res.multisplit <- multisplit(x = sim.geno, y = y, family = "gaussian", B = B)

  # test hierarchy: Tippett
  set.seed(2)
  res.T <- test_hierarchy(x = sim.geno, y = y, dendr = dendr, family = "gaussian",
                          B = B)

  # test hierarchy: Stouffer
  set.seed(2)
  res.S <- test_hierarchy(x = sim.geno, y = y, dendr = dendr, family = "gaussian",
                          B = B, agg.method = "Stouffer")

  ## Test
  # This list encodes the tree structure
  cluster_test <- list(c("rsid1", "rsid2", "rsid3", "rsid4", "rsid5"),
                       c("rsid1", "rsid2"),
                       c("rsid3", "rsid4", "rsid5"),
                       c("rsid3", "rsid4"),
                       "rsid1",
                       "rsid2",
                       "rsid3",
                       "rsid4",
                       "rsid5")

  pvals_to_be <- check_test_hierarchy(x = sim.geno, y = y, clvar = NULL,
                                      res.multisplit = res.multisplit,
                                      B = B, cluster_test = cluster_test)



  # Tippett
  compare_with <- pvals_to_be

  expected_result <- data.frame(block = c(NA, NA),
                                p.value = compare_with[c("rsid1", "rsid4")])
  expected_result$significant.cluster <- list(c("rsid1"), c("rsid4"))
  rownames(expected_result) <- NULL
  attr(expected_result, "class") <- c("data.frame")

  expect_equal(res.T$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-120)
  expect_equal(res.T$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-90)
  expect_equal(res.T$res.hierarchy$significant.cluster,
               expected_result$significant.cluster)

  # Stouffer
  compare_with <- pvals_to_be

  expected_result <- data.frame(block = c(NA, NA),
                                p.value = compare_with[c("rsid1", "rsid4")])
  expected_result$significant.cluster <- list(c("rsid1"), c("rsid4"))
  rownames(expected_result) <- NULL
  attr(expected_result, "class") <- c("data.frame")

  expect_equal(res.S$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-120)
  expect_equal(res.S$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-90)
  expect_equal(res.S$res.hierarchy$significant.cluster,
               expected_result$significant.cluster)
})

### Example II ###
test_that("test_hierarchy: check output (Example II)", {
  ## simulate index
  n <- 800
  p <- 5
  B <- 5

  ## simulate data
  require(MASS)
  set.seed(9229)
  sim.geno <- mvrnorm(n = n, mu = rep(0, p), Sigma = toeplitz(0.8^(seq(0, p - 1))))
  colnames(sim.geno) <- paste0("rsid", 1:p)
  sim.clvar <- matrix(rnorm(n * 3), ncol = 3)
  colnames(sim.clvar) <- paste0("clvar", 1:3)

  set.seed(144)
  data.dim <- dim(sim.geno) # first entry corresponds to rows and second to columns

  ind.active <- sample(1:data.dim[2], 2)
  beta <- rep(0, data.dim[2])
  beta[ind.active] <- 2
  y <- sim.geno %*% beta + sim.clvar %*% c(0.25, 0.5, 1) + rnorm(data.dim[1])

  # cluster the data
  dendr <- cluster_var(x = sim.geno)
  # plot(dendr$res.tree[[1]])

  # multisplit
  set.seed(555)
  res.multisplit <- multisplit(x = sim.geno, y = y, clvar = sim.clvar,
                               family = "gaussian", B = B)
  # test hierarchy: Tippett
  set.seed(555)
  res.T <- test_hierarchy(x = sim.geno, y = y, clvar = sim.clvar,
                          dendr = dendr, family = "gaussian", B = B)

  # test hierarchy: Stouffer
  set.seed(555)
  res.S <- test_hierarchy(x = sim.geno, y = y, clvar = sim.clvar,
                          dendr = dendr, family = "gaussian", B = B,
                          agg.method = "Stouffer")

  ## test
  # This list encodes the tree structure
  cluster_test <- list(c("rsid1", "rsid2", "rsid3", "rsid4", "rsid5"),
                       c("rsid1", "rsid2"),
                       c("rsid3", "rsid4", "rsid5"),
                       c("rsid3", "rsid4"),
                       "rsid1",
                       "rsid2",
                       "rsid3",
                       "rsid4",
                       "rsid5")

  pvals_to_be <- check_test_hierarchy(x = sim.geno, y = y, clvar = sim.clvar,
                                      res.multisplit = res.multisplit,
                                      B = B, cluster_test = cluster_test)

  # Tippett
  compare_with <- pvals_to_be

  expected_result <- data.frame(block = c(NA, NA),
                                p.value = compare_with[c("rsid1", "rsid4")])
  expected_result$significant.cluster <- list(c("rsid1"), c("rsid4"))
  rownames(expected_result) <- NULL
  attr(expected_result, "class") <- c("data.frame")

  expect_equal(res.T$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-120)
  expect_equal(res.T$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-90)
  expect_equal(res.T$res.hierarchy$significant.cluster,
               expected_result$significant.cluster)

  # Stouffer
  compare_with <- pvals_to_be

  expected_result <- data.frame(block = c(NA, NA),
                                p.value = compare_with[c("rsid1", "rsid4")])
  expected_result$significant.cluster <- list(c("rsid1"), c("rsid4"))
  rownames(expected_result) <- NULL
  attr(expected_result, "class") <- c("data.frame")

  expect_equal(res.S$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-115)
  expect_equal(res.S$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-85)
  expect_equal(res.S$res.hierarchy$significant.cluster,
               expected_result$significant.cluster)
})


#### Check output with multiple data sets ####
# Function for generating the data
require(MASS)

gen_one <- function(n, p, seed1, seed2, seed3, num_clvar = NULL,
                    coef_clvar = NULL) {
  set.seed(seed1)
  x <- mvrnorm(n = n, mu = rep(0, p), Sigma = toeplitz(0.8^(seq(0, p - 1))) )
  colnames(x) <- paste0("rsid", 1:p)
  if (!is.null(num_clvar)) {
    clvar <- matrix(rnorm(n * 3), ncol = 3)
    colnames(clvar) <- paste0("clvar", 1:3)
  } else {
    clvar <- NULL
  }

  set.seed(seed2)
  data.dim <- dim(x) # first entry corresponds to rows and second to columns
  ind.active <- sample(1:data.dim[2], 2)
  beta <- rep(0, data.dim[2])
  beta[ind.active] <- 2

  set.seed(seed3)
  if (!is.null(num_clvar)) {
    y <- x %*% beta + rnorm(data.dim[1]) + clvar %*% coef_clvar
  } else {
    y <- x %*% beta + rnorm(data.dim[1])
  }


  return(list(x = x, y = y, clvar = clvar))
}

### Example III ###
test_that("test_hierarchy: check output (Example III multiple data sets)", {
  skip_on_bioc()

  ## simulate index
  n <- 800
  p <- 5
  B <- 50

  ## simulate data
  r1 <- gen_one(n = n, p = p, seed1 = 9229, seed2 = 144, seed3 = 8)
  r2 <- gen_one(n = n, p = p, seed1 = 929, seed2 = 144, seed3 = 99)
  r3 <- gen_one(n = n, p = p, seed1 = 99, seed2 = 144, seed3 = 100)
  r4 <- gen_one(n = n, p = p, seed1 = 9, seed2 = 144, seed3 = 1111)

  x <- list(r1$x, r2$x, r3$x, r4$x)
  y <- list(r1$y, r2$y, r3$y, r4$y)
  # clvar <- list(r1$clvar, r2$clvar, r3$clvar, r4$clvar)

  # cluster the data
  dendr <- cluster_var(x = x)
  # plot(dendr$res.tree[[1]])

  # multisplit
  set.seed(744)
  res.multisplit <- multisplit(x = x, y = y, family = "gaussian", B = B)

  # test hierarchy: Tippett
  set.seed(744)
  res.T <- test_hierarchy(x = x, y = y, dendr = dendr,
                          family = "gaussian", B = B)

  # test hierarchy: Stouffer
  set.seed(744)
  res.S <- test_hierarchy(x = x, y = y, dendr = dendr,
                          family = "gaussian", B = B,
                          agg.method = "Stouffer")

  ## Test
  # This list encodes the tree structure
  cluster_test <- list(c("rsid1", "rsid2", "rsid3", "rsid4", "rsid5"),
                       c("rsid1", "rsid2", "rsid3"),
                       c("rsid4", "rsid5"),
                       c("rsid1", "rsid2"),
                       "rsid1",
                       "rsid2",
                       "rsid3",
                       "rsid4",
                       "rsid5")

  res1 <- check_test_hierarchy(x = x[[1]], y = y[[1]], clvar = NULL,
                               res.multisplit = res.multisplit[1],
                               B = B, cluster_test = cluster_test)

  res2 <- check_test_hierarchy(x = x[[2]], y = y[[2]], clvar = NULL,
                               res.multisplit = res.multisplit[2],
                               B = B, cluster_test = cluster_test)

  res3 <- check_test_hierarchy(x = x[[3]], y = y[[3]], clvar = NULL,
                               res.multisplit = res.multisplit[3],
                               B = B, cluster_test = cluster_test)

  res4 <- check_test_hierarchy(x = x[[4]], y = y[[4]], clvar = NULL,
                               res.multisplit = res.multisplit[4],
                               B = B, cluster_test = cluster_test)

  pvals_to_be <- rbind(res1, res2, res3, res4)

  # Tippett
  compare_with <- apply(X = pvals_to_be, MARGIN = 2,
                        FUN = function(x, len_y) {
                          max(1 - (1 - min(x))^(len_y), .Machine$double.neg.eps)
                        },
                        len_y = 4)

  expected_result <- data.frame(block = c(NA, NA),
                                p.value = compare_with[c("rsid4", "rsid1")])
  expected_result$significant.cluster <- list(c("rsid4"), c("rsid1"))
  rownames(expected_result) <- NULL
  attr(expected_result, "class") <- c("data.frame")

  expect_equal(res.T$res.hierarchy$p.value, expected_result$p.value,
               tol = 1e-145)
  expect_equal(res.T$res.hierarchy$significant.cluster,
               expected_result$significant.cluster)

  # Stouffer
  # c(0.5, 0.5, 0.5, 0.5)
  stouffer_weights <- sqrt(c(800, 800, 800, 800) / sum(c(800, 800, 800, 800)))

  compare_with <- apply(X = pvals_to_be, MARGIN = 2,
                        FUN = function(x, len_y, stouffer_weights) {
                          pnorm(sum(stouffer_weights * qnorm(x)))
                        },
                        len_y = 4, stouffer_weights = stouffer_weights)


  expected_result <- data.frame(block = c(NA, NA),
                                p.value = compare_with[c("rsid4", "rsid1")])
  expected_result$significant.cluster <- list(c("rsid4"), c("rsid1"))
  rownames(expected_result) <- NULL
  attr(expected_result, "class") <- c("data.frame")

  expect_equal(res.S$res.hierarchy$p.value, expected_result$p.value,
               tol = 1e-200)
  expect_equal(res.S$res.hierarchy$significant.cluster,
               expected_result$significant.cluster)
})

### Example IV unbalanced data sets ###
test_that("test_hierarchy: check output (Example IV multiple data sets)", {
  skip_on_bioc()

  ## simulate index
  p <- 5
  B <- 50

  ## simulate data
  require(MASS)

  r1 <- gen_one(n = 800, p = p, seed1 = 9229, seed2 = 144, seed3 = 8)
  r2 <- gen_one(n = 200, p = p, seed1 = 929, seed2 = 144, seed3 = 99)
  r3 <- gen_one(n = 350, p = p, seed1 = 99, seed2 = 144, seed3 = 100)
  r4 <- gen_one(n = 50, p = p, seed1 = 9, seed2 = 144, seed3 = 1111)

  x <- list(r1$x, r2$x, r3$x, r4$x)
  y <- list(r1$y, r2$y, r3$y, r4$y)
  # clvar <- list(r1$clvar, r2$clvar, r3$clvar, r4$clvar)

  # cluster the data
  dendr <- cluster_var(x = x)
  # plot(dendr$res.tree[[1]])

  # multisplit
  set.seed(6)
  res.multisplit <- multisplit(x = x, y = y, family = "gaussian", B = B)

  # test hierarchy: Tippett
  set.seed(6)
  res.T <- test_hierarchy(x = x, y = y, dendr = dendr,
                          family = "gaussian", B = B)

  # test hierarchy: Stouffer
  set.seed(6)
  res.S <- test_hierarchy(x = x, y = y, dendr = dendr,
                          family = "gaussian", B = B,
                          agg.method = "Stouffer")


  ## Test
  # This list encodes the tree structure
  cluster_test <- list(c("rsid1", "rsid2", "rsid3", "rsid4", "rsid5"),
                       c("rsid1", "rsid2", "rsid3"),
                       c("rsid4", "rsid5"),
                       c("rsid1", "rsid2"),
                       "rsid1",
                       "rsid2",
                       "rsid3",
                       "rsid4",
                       "rsid5")

  res1 <- check_test_hierarchy(x = x[[1]], y = y[[1]], clvar = NULL,
                               res.multisplit = res.multisplit[1],
                               B = B, cluster_test = cluster_test)

  res2 <- check_test_hierarchy(x = x[[2]], y = y[[2]], clvar = NULL,
                               res.multisplit = res.multisplit[2],
                               B = B, cluster_test = cluster_test)

  res3 <- check_test_hierarchy(x = x[[3]], y = y[[3]], clvar = NULL,
                               res.multisplit = res.multisplit[3],
                               B = B, cluster_test = cluster_test)

  res4 <- check_test_hierarchy(x = x[[4]], y = y[[4]], clvar = NULL,
                               res.multisplit = res.multisplit[4],
                               B = B, cluster_test = cluster_test)

  pvals_to_be <- rbind(res1, res2, res3, res4)

  # Tippett
  compare_with <- apply(X = pvals_to_be, MARGIN = 2,
                        FUN = function(x, len_y) {
                          max(1 - (1 - min(x))^(len_y), .Machine$double.neg.eps)
                        },
                        len_y = 4)

  expected_result <- data.frame(block = c(NA, NA),
                                p.value = compare_with[c("rsid1", "rsid4")])
  expected_result$significant.cluster <- list(c("rsid1"), c("rsid4"))
  rownames(expected_result) <- NULL

  expect_equal(res.T$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-140)
  expect_equal(res.T$res.hierarchy$p.value[2], expected_result$p.value[2],
               tol = 1e-140)
  expect_equal(res.T$res.hierarchy$significant.cluster,
               expected_result$significant.cluster)

  # Stouffer
  stouffer_weights <- sqrt(c(800, 200, 350, 50) / sum(c(800, 200, 350, 50)))

  compare_with <- apply(X = pvals_to_be, MARGIN = 2,
                        FUN = function(x, len_y, stouffer_weights) {
                          pnorm(sum(stouffer_weights * qnorm(x)))
                        },
                        len_y = 4, stouffer_weights = stouffer_weights)

  expected_result <- data.frame(block = c(NA, NA),
                                p.value = compare_with[c("rsid1", "rsid4")])
  expected_result$significant.cluster <- list(c("rsid1"), c("rsid4"))
  rownames(expected_result) <- NULL

  expect_equal(res.S$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-220)
  expect_equal(res.S$res.hierarchy$p.value[2], expected_result$p.value[2],
               tol = 1e-190)
  expect_equal(res.S$res.hierarchy$significant.cluster,
               expected_result$significant.cluster)
})

### Example V with co-variables ###
test_that("test_hierarchy: check output (Example V multiple data sets)", {
  skip_on_bioc()

  ## simulate index
  p <- 5
  B <- 50

  ## simulate data
  require(MASS)

  r1 <- gen_one(n = 800, p = p, seed1 = 9229, seed2 = 144, seed3 = 8,
                num_clvar = 3, coef_clvar = c(0.5, 0.25, 1.25))
  r2 <- gen_one(n = 200, p = p, seed1 = 929, seed2 = 144, seed3 = 99,
                num_clvar = 3, coef_clvar = c(0.5, 0.25, 1.25))
  r3 <- gen_one(n = 350, p = p, seed1 = 99, seed2 = 144, seed3 = 100,
                num_clvar = 3, coef_clvar = c(0.5, 0.25, 1.25))
  r4 <- gen_one(n = 50, p = p, seed1 = 9, seed2 = 144, seed3 = 1111,
                num_clvar = 3, coef_clvar = c(0.5, 0.25, 1.25))

  x <- list(r1$x, r2$x, r3$x, r4$x)
  y <- list(r1$y, r2$y, r3$y, r4$y)
  clvar <- list(r1$clvar, r2$clvar, r3$clvar, r4$clvar) # with co-variables

  # cluster the data
  dendr <- cluster_var(x = x)
  # plot(dendr$res.tree[[1]])

  # multisplit
  set.seed(3)
  res.multisplit <- multisplit(x = x, y = y, clvar = clvar, family = "gaussian",
                               B = B)

  # test hierarchy: Tippett
  set.seed(3)
  res.T <- test_hierarchy(x = x, y = y, clvar = clvar,
                          dendr = dendr, family = "gaussian",
                          B = B)

  # test hierarchy: Stouffer
  set.seed(3)
  res.S <- test_hierarchy(x = x, y = y, clvar = clvar,
                          dendr = dendr, family = "gaussian",
                          B = B, agg.method = "Stouffer")

  ## Test
  # This list encodes the tree structure
  cluster_test <- list(c("rsid1", "rsid2", "rsid3", "rsid4", "rsid5"),
                       c("rsid1", "rsid2"),
                       c("rsid3", "rsid4", "rsid5"),
                       c("rsid4", "rsid5"),
                       "rsid1",
                       "rsid2",
                       "rsid3",
                       "rsid4",
                       "rsid5")

  res1 <- check_test_hierarchy(x = x[[1]], y = y[[1]], clvar = clvar[[1]],
                               res.multisplit = res.multisplit[1],
                               B = B, cluster_test = cluster_test)

  res2 <- check_test_hierarchy(x = x[[2]], y = y[[2]], clvar = clvar[[2]],
                               res.multisplit = res.multisplit[2],
                               B = B, cluster_test = cluster_test)

  res3 <- check_test_hierarchy(x = x[[3]], y = y[[3]], clvar = clvar[[3]],
                               res.multisplit = res.multisplit[3],
                               B = B, cluster_test = cluster_test)

  res4 <- check_test_hierarchy(x = x[[4]], y = y[[4]], clvar = clvar[[4]],
                               res.multisplit = res.multisplit[4],
                               B = B, cluster_test = cluster_test)

  pvals_to_be <- rbind(res1, res2, res3, res4)

  # Tippett
  compare_with <- apply(X = pvals_to_be, MARGIN = 2,
                        FUN = function(x, len_y) {
                          max(1 - (1 - min(x))^(len_y), .Machine$double.neg.eps)
                        },
                        len_y = 4)

  expected_result <- data.frame(block = c(NA, NA),
                                p.value = compare_with[c("rsid1", "rsid4")])
  expected_result$significant.cluster <- list(c("rsid1"), c("rsid4"))
  rownames(expected_result) <- NULL
  attr(expected_result, "class") <- c("data.frame")

  expect_equal(res.T$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-140)
  expect_equal(res.T$res.hierarchy$p.value[2], expected_result$p.value[2],
               tol = 1e-140)
  expect_equal(res.T$res.hierarchy$significant.cluster,
               expected_result$significant.cluster)

  # Stouffer
  stouffer_weights <- sqrt(c(800, 200, 350, 50) / sum(c(800, 200, 350, 50)))

  compare_with <- apply(X = pvals_to_be, MARGIN = 2,
                        FUN = function(x, len_y, stouffer_weights) {
                          pnorm(sum(stouffer_weights * qnorm(x)))
                        },
                        len_y = 4, stouffer_weights = stouffer_weights)

  expected_result <- data.frame(block = c(NA, NA),
                                p.value = compare_with[c("rsid1", "rsid4")])
  expected_result$significant.cluster <- list(c("rsid1"), c("rsid4"))
  rownames(expected_result) <- NULL
  attr(expected_result, "class") <- c("data.frame")

  expect_equal(res.S$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-220)
  expect_equal(res.S$res.hierarchy$p.value[2], expected_result$p.value[2],
               tol = 1e-180)
  expect_equal(res.S$res.hierarchy$significant.cluster,
               expected_result$significant.cluster)


  ### Example VI, same data as in Example V but without clvar ###
  # test hierarchy: no clvar but the response was created incl. clvar
  # test_hierarchy: check output (Example VI multiple data sets)

  # multisplit
  set.seed(4)
  res.multisplit <- multisplit(x = x, y = y, family = "gaussian",
                               B = B)

  # test hierarchy: Tippett
  set.seed(4)
  res.T <- test_hierarchy(x = x, y = y, dendr = dendr,
                          family = "gaussian", B = B)

  # test hierarchy: Stouffer
  set.seed(4)
  res.S <- test_hierarchy(x = x, y = y, dendr = dendr,
                          family = "gaussian", B = B,
                          agg.method = "Stouffer")

  cluster_test <- list(c("rsid1", "rsid2", "rsid3", "rsid4", "rsid5"),
                       c("rsid1", "rsid2"),
                       c("rsid3", "rsid4", "rsid5"),
                       c("rsid4", "rsid5"),
                       "rsid1",
                       "rsid2",
                       "rsid3",
                       "rsid4",
                       "rsid5")

  res1 <- check_test_hierarchy(x = x[[1]], y = y[[1]], clvar = NULL,
                               res.multisplit = res.multisplit[1],
                               B = B, cluster_test = cluster_test)

  res2 <- check_test_hierarchy(x = x[[2]], y = y[[2]], clvar = NULL,
                               res.multisplit = res.multisplit[2],
                               B = B, cluster_test = cluster_test)

  res3 <- check_test_hierarchy(x = x[[3]], y = y[[3]], clvar = NULL,
                               res.multisplit = res.multisplit[3],
                               B = B, cluster_test = cluster_test)

  res4 <- check_test_hierarchy(x = x[[4]], y = y[[4]], clvar = NULL,
                               res.multisplit = res.multisplit[4],
                               B = B, cluster_test = cluster_test)

  pvals_to_be <- rbind(res1, res2, res3, res4)

  # Tippett
  compare_with <- apply(X = pvals_to_be, MARGIN = 2,
                        FUN = function(x, len_y) {
                          max(1 - (1 - min(x))^(len_y), .Machine$double.neg.eps)
                        },
                        len_y = 4)

  expected_result <- data.frame(block = c(NA, NA),
                                p.value = compare_with[c("rsid1", "rsid4")])
  expected_result$significant.cluster <- list(c("rsid1"), c("rsid4"))
  rownames(expected_result) <- NULL

  expect_equal(res.T$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-100)
  expect_equal(res.T$res.hierarchy$p.value[2], expected_result$p.value[2],
               tol = 1e-100)
  expect_equal(res.T$res.hierarchy$significant.cluster,
               expected_result$significant.cluster)

  # Stouffer
  stouffer_weights <- sqrt(c(800, 200, 350, 50) / sum(c(800, 200, 350, 50)))

  compare_with <- apply(X = pvals_to_be, MARGIN = 2,
                        FUN = function(x, len_y, stouffer_weights) {
                          pnorm(sum(stouffer_weights * qnorm(x)))
                        },
                        len_y = 4, stouffer_weights = stouffer_weights)

  expected_result <- data.frame(block = c(NA, NA),
                                p.value = compare_with[c("rsid3_rsid4_rsid5", "rsid1")])
  expected_result$significant.cluster <- list(c("rsid3", "rsid4", "rsid5"), c("rsid1"))
  rownames(expected_result) <- NULL

  expect_equal(res.S$res.hierarchy$p.value[1], expected_result$p.value[1],
               tol = 1e-112)
  expect_equal(res.S$res.hierarchy$p.value[2], expected_result$p.value[2],
               tol = 1e-100)
  expect_equal(res.S$res.hierarchy$significant.cluster,
               expected_result$significant.cluster)
})

