require("testthat")

test_that("multisplit: check input", {
  expect_error(multisplit(x = NULL, y = NULL, family = "gaussian"),
               "The response y is required to be a vector or a list of vectors if multiple data sets are present.")

  # Please note that family = NULL results in taking the default value.
  expect_error(multisplit(x = NULL, y = matrix(1:4, ncol = 2), family = NULL),
               "The elements of the list y are required to be numeric vectors or matrices with only one column. In the case of only one data set, it is enough that y is a numeric vector or matrix with only one column but it can as well be a list with one element.")

  expect_error(multisplit(x = NULL, y = list(matrix(1:4, ncol = 2), matrix(1:4, ncol = 2)), family = NULL),
               "The elements of the list y are required to be numeric vectors or matrices with only one column.")

  expect_error(multisplit(x = NULL, y = 1:2, family = "gaussian"),
               "The input x is required to be a matrix or a list of matrices if multiple data sets are present.")

  expect_error(multisplit(x = matrix(1:4, ncol = 2), y = 1:2, family = "gaussian"),
               "The matrix x is required to have column names. If there is no natural naming convention, then one can set them to some integer, say, 1 to p.",
               fixed = TRUE)

  set.seed(88)
  x <- matrix(rnorm(100), ncol = 2)
  colnames(x) <- c("col1", "col2")
  y <- 1:50
  # There occures an error because there are not enough observations for glmnet
  expect_warning(multisplit(x = x[1:3, ], y = y[1:3], family = "gaussian",
                            B = 20),
                 "There occurred some errors while multi-sample splitting. See attribute 'errorMsgs' of the return object for more details.")
  # TODO maybe return the warnings and errors combined for each data set

  expect_error(multisplit(x = x, y = y, family = "gaussian", B = NULL),
               "The argument B has to be a positive integer.")

  expect_error(multisplit(x = x, y = y, family = "gaussian", B = -2),
               "The argument B has to be a positive integer.")

})



check_dim <- function(n, p, B, proportion.select){
  sim.geno <- matrix(rnorm(n * p), ncol = p, nrow = n)
  colnames(sim.geno) <- paste0("lll", 1:p)
  sim.coef <- rnorm(p, sd = 0.25)
  y <- sim.geno %*% sim.coef + rnorm(10, sd = 0.5)
  res_multisplit <- multisplit(x = sim.geno, y = y, B = B,
                               proportion.select = proportion.select,
                               family = "gaussian")
  return(list(x = sim.geno, res_multisplit = res_multisplit))
}

test_that("multisplit: check return object", {
  # Check if the dimensions of the output of multisplit is correct.
  n <- 10; p <- 10; B <- 20; proportion.select <- 1/6
  res_tmp <- check_dim(n = n, p = p, B = B,
                              proportion.select = proportion.select)
  res_multisplit <- res_tmp$res_multisplit; x <- res_tmp$x
  expect_true(all(dim(res_multisplit[[1]]$out.sample) == c(B, n / 2)))
  expect_true(all(dim(res_multisplit[[1]]$sel.coef) ==
                c(B, floor(n * proportion.select))))
  sel.coef <- res_multisplit[[1]]$sel.coef
  expect_true(all(sel.coef[!is.na(sel.coef)] %in% colnames(x)))

  n <- 1000; p <- 100; B <- 20; proportion.select <- 1/6
  res_tmp <- check_dim(n = n, p = p, B = B,
                       proportion.select = proportion.select)
  res_multisplit <- res_tmp$res_multisplit; x <- res_tmp$x
  expect_true(all(dim(res_multisplit[[1]]$out.sample) == c(B, n / 2)))
  expect_true(all(dim(res_multisplit[[1]]$sel.coef) ==
                    c(B, floor(n * proportion.select))))
  sel.coef <- res_multisplit[[1]]$sel.coef
  expect_true(all(sel.coef[!is.na(sel.coef)] %in% colnames(x)))
  # If there are less variables than number of selected variables, then there
  # have to be NA's in the matrix of selected variables.
  if ((tmp <- p - floor(n * proportion.select)) < 0) {
    expect_true(all(rowSums(is.na(sel.coef)) >= abs(tmp)))
  }

  n <- 100; p <- 1000; B <- 20; proportion.select <- 1/6
  res_tmp <- check_dim(n = n, p = p, B = B,
                       proportion.select = proportion.select)
  res_multisplit <- res_tmp$res_multisplit; x <- res_tmp$x
  expect_true(all(dim(res_multisplit[[1]]$out.sample) == c(B, n / 2)))
  expect_true(all(dim(res_multisplit[[1]]$sel.coef) ==
                    c(B, floor(n * proportion.select))))
  sel.coef <- res_multisplit[[1]]$sel.coef
  expect_true(all(sel.coef[!is.na(sel.coef)] %in% colnames(x)))
})

## TODO Run this test only locally. (Not suitable for Windows)

# test_that("multisplit: check if the functions runs in parallel", {
#   n <- 100; p <- 1000; B <- 20; proportion.select <- 1 / 6
#   sim.geno1 <- matrix(rnorm(n * p), ncol = p, nrow = n)
#   colnames(sim.geno1) <- paste0("lll", 1:p)
#
#   sim.geno2 <- matrix(rnorm(n * p), ncol = p, nrow = n)
#   colnames(sim.geno2) <- paste0("lll", 1:p)
#
#   sim.coef <- rnorm(p, sd = 0.25)
#   y1 <- sim.geno1 %*% sim.coef + rnorm(10, sd = 0.5)
#   y2 <- sim.geno2 %*% sim.coef + rnorm(10, sd = 0.5)
#   res_multisplit <- multisplit(x = list(sim.geno1, sim.geno2), y = list(y1, y2),
#                                proportion.select = proportion.select, B = B,
#                                family = "gaussian", seed = 789,
#                                parallel = "multicore", ncpus = 2)
#   expect_true(length(res_multisplit) == 2)
#   expect_is(res_multisplit, "hierM")
#   expect_true(all(dim(res_multisplit[[1]]$out.sample) == c(B, n / 2)))
#   expect_true(all(dim(res_multisplit[[2]]$out.sample) == c(B, n / 2)))
#   expect_true(all(dim(res_multisplit[[1]]$sel.coef) ==
#                     c(B, floor(n * proportion.select))))
#   expect_true(all(dim(res_multisplit[[2]]$sel.coef) ==
#                     c(B, floor(n * proportion.select))))
#   sel.coef <- res_multisplit[[1]]$sel.coef
#   expect_true(all(sel.coef[!is.na(sel.coef)] %in% colnames(sim.geno1)))
#   sel.coef <- res_multisplit[[2]]$sel.coef
#   expect_true(all(sel.coef[!is.na(sel.coef)] %in% colnames(sim.geno2)))
# })

