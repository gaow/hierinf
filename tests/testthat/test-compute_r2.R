require("testthat")

## random number generator
RNGkind("L'Ecuyer-CMRG")

test_that("compute_r2: check input", {
  expect_error(compute_r2(x = NULL, y = NULL, res.test.hierarchy = NULL,
                          family = NULL),
               "The argument res.test.hierarchy is required to be the output of the function test_hierarchy, the function test_only_hierarchy, or the function multisplit.")

  set.seed(88)
  x <- matrix(rnorm(100), ncol = 2)
  colnames(x) <- c("col1", "col2")
  y <- 1:50
  res.multisplit <- multisplit(x = x, y = y, B = 5)
  expect_error(compute_r2(x = NULL, y = NULL,
                                   res.test.hierarchy = res.multisplit,
                                   family = NULL),
               "The response y is required to be a vector or a list of vectors if multiple data sets are present.")

  expect_error(compute_r2(x = NULL, y = y,
                          res.test.hierarchy = res.multisplit,
                          family = NULL),
               "The input x is required to be a matrix or a list of matrices if multiple data sets are present.")

  x <- 1:50
  expect_error(compute_r2(x = x, y = y,
                          res.test.hierarchy = res.multisplit,
                          family = NULL),
               "The input x is required to be a matrix or a list of matrices if multiple data sets are present.")

  x <- matrix(1:50, ncol = 1)
  expect_error(compute_r2(x = x, y = y,
                          res.test.hierarchy = res.multisplit,
                          family = NULL),
               "The matrix x is required to have column names. If there is no natural naming convention, then one can set them to some integer, say, 1 to p.")

  x <- matrix(1:50, ncol = 1)
  colnames(x) <- "CCOOLL1"
  expect_error(compute_r2(x = x, y = y,
                          res.test.hierarchy = res.multisplit,
                          family = NULL),
               "The selected variables in the output of the function call to test_hierarchy or multisplit does not match the column names of the argument x.")

  set.seed(88)
  x <- matrix(rnorm(100), ncol = 2)
  colnames(x) <- c("col1", "col2")
  y <- 1:50
  res.multisplit <- multisplit(x = x, y = y, B = 5)
  expect_error(compute_r2(x = x, y = y,
                          res.test.hierarchy = res.multisplit,
                          family = NULL, colnames.cluster = "CC_col"),
               "Each variable which column names is specified in the argument colnames.cluster has to be contained in the data set or at least one data set for multiple data sets.")

 # Please note that family = NULL results in takeing the default value.
})

#### Check output with simulated data ####
test_that("compute_r2: check output", {
  ### Example I noise ###
  set.seed(88)
  x <- matrix(rnorm(100), ncol = 2)
  colnames(x) <- c("col1", "col2")
  y <- 1:50
  res.multisplit <- multisplit(x = x, y = y, B = 5)

  dat <- data.frame(y = y, x)
  out <- res.multisplit[[1]]$out.sample
  sel <- res.multisplit[[1]]$sel.coef

  out1 <- out[1, ]
  sel1 <- sel[1, ][!is.na(sel[1, ])]
  r1 <- summary(lm(y ~ ., data = dat[out1, c("y", sel1)]))$adj.r.squared

  out2 <- out[2, ]
  sel2 <- sel[2, ][!is.na(sel[2, ])]
  r2 <- summary(lm(y ~ ., data = dat[out2, c("y", sel2)]))$adj.r.squared

  out3 <- out[3, ]
  sel3 <- sel[3, ][!is.na(sel[3, ])]
  r3 <- summary(lm(y ~ ., data = dat[out3, c("y", sel3)]))$adj.r.squared

  out4 <- out[4, ]
  sel4 <- sel[4, ][!is.na(sel[4, ])]
  r4 <- summary(lm(y ~ ., data = dat[out4, c("y", sel4)]))$adj.r.squared

  out5 <- out[5, ]
  sel5 <- sel[5, ][!is.na(sel[5, ])]
  r5 <- summary(lm(y ~ ., data = dat[out5, c("y", sel5)]))$adj.r.squared

  expect_equal(compute_r2(x = x, y = y,
                          res.test.hierarchy = res.multisplit,
                          family = NULL), mean(c(r1, r2, r3, r4, r5)))


  ### Example II with signal ###
  set.seed(88)
  x <- matrix(rnorm(100), ncol = 2)
  colnames(x) <- c("col1", "col2")
  y <- x %*% c(2, 2) + rnorm(50)
  res.multisplit <- multisplit(x = x, y = y, B = 5)

  dat <- data.frame(y = y, x)
  out <- res.multisplit[[1]]$out.sample
  sel <- res.multisplit[[1]]$sel.coef

  out1 <- out[1, ]
  sel1 <- sel[1, ][!is.na(sel[1, ])]
  r1 <- summary(lm(y ~ ., data = dat[out1, c("y", sel1)]))$adj.r.squared

  out2 <- out[2, ]
  sel2 <- sel[2, ][!is.na(sel[2, ])]
  r2 <- summary(lm(y ~ ., data = dat[out2, c("y", sel2)]))$adj.r.squared

  out3 <- out[3, ]
  sel3 <- sel[3, ][!is.na(sel[3, ])]
  r3 <- summary(lm(y ~ ., data = dat[out3, c("y", sel3)]))$adj.r.squared

  out4 <- out[4, ]
  sel4 <- sel[4, ][!is.na(sel[4, ])]
  r4 <- summary(lm(y ~ ., data = dat[out4, c("y", sel4)]))$adj.r.squared

  out5 <- out[5, ]
  sel5 <- sel[5, ][!is.na(sel[5, ])]
  r5 <- summary(lm(y ~ ., data = dat[out5, c("y", sel5)]))$adj.r.squared

  expect_equal(compute_r2(x = x, y = y,
                          res.test.hierarchy = res.multisplit,
                          family = NULL), mean(c(r1, r2, r3, r4, r5)))

  ### Example III multiple data sets (unbalanced) ###
  set.seed(88)
  # n1 = 50
  x1 <- matrix(rnorm(100), ncol = 2)
  colnames(x1) <- c("col1", "col2")
  y1 <- x1 %*% c(2, 2) + rnorm(50)

  # n2 = 50
  x2 <- matrix(rnorm(200), ncol = 2)
  colnames(x2) <- c("col1", "col2")
  y2 <- x2 %*% c(4, 4) + rnorm(100)

  res.multisplit <- multisplit(x = list(x1, x2), y = list(y1, y2), B = 3)

  dat1 <- data.frame(y = y1, x1)
  out <- res.multisplit[[1]]$out.sample
  sel <- res.multisplit[[1]]$sel.coef

  out1 <- out[1, ]
  sel1 <- sel[1, ][!is.na(sel[1, ])]
  r1.1 <- summary(lm(y ~ ., data = dat1[out1, c("y", sel1)]))$adj.r.squared

  out2 <- out[2, ]
  sel2 <- sel[2, ][!is.na(sel[2, ])]
  r1.2 <- summary(lm(y ~ ., data = dat1[out2, c("y", sel2)]))$adj.r.squared

  out3 <- out[3, ]
  sel3 <- sel[3, ][!is.na(sel[3, ])]
  r1.3 <- summary(lm(y ~ ., data = dat1[out3, c("y", sel3)]))$adj.r.squared


  dat2 <- data.frame(y = y2, x2)
  out <- res.multisplit[[2]]$out.sample
  sel <- res.multisplit[[2]]$sel.coef

  out1 <- out[1, ]
  sel1 <- sel[1, ][!is.na(sel[1, ])]
  r2.1 <- summary(lm(y ~ ., data = dat2[out1, c("y", sel1)]))$adj.r.squared

  out2 <- out[2, ]
  sel2 <- sel[2, ][!is.na(sel[2, ])]
  r2.2 <- summary(lm(y ~ ., data = dat2[out2, c("y", sel2)]))$adj.r.squared

  out3 <- out[3, ]
  sel3 <- sel[3, ][!is.na(sel[3, ])]
  r2.3 <- summary(lm(y ~ ., data = dat2[out3, c("y", sel3)]))$adj.r.squared


  expect_equal(compute_r2(x = list(x1, x2), y = list(y1, y2),
                          res.test.hierarchy = res.multisplit,
                          family = NULL),
               100 / 300 * mean(c(r1.1, r1.2, r1.3)) +
                 200 / 300 * mean(c(r2.1, r2.2, r2.3)))
})

