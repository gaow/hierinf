#' Compute R squared
#'
#' Compute the R squared value for a given cluster or group of variables.
#'
#' @param x a matrix or list of matrices for multiple data sets. The matrix or
#' matrices have to be of type numeric and are required to have column names
#' / variable names. The rows and the columns represent the observations and
#' the variables, respectively.
#' @param y a vector, a matrix with one column, or list of the aforementioned
#' objects for multiple data sets. The vector, vectors, matrix, or matrices
#' have to be of type numeric.
#' @param res.test.hierarchy the output of one of the functions
#' \code{\link[hierinf]{test_hierarchy}},
#' \code{\link[hierinf]{test_only_hierarchy}}, or
#' \code{\link[hierinf]{multisplit}}.
#' @param clvar a matrix or list of matrices of control variables.
#' @param family a character string naming a family of the error distribution;
#' either \code{"gaussian"} or \code{"binomial"}.
#' @param colnames.cluster The column names / variables names of the cluster
#' of interest. If not supplied, the R squared value of the full model is
#' computed.
#'
#' @details The R squared value is computed based on the output of the multi-sample
#' splitting step. For each split, the intersection of the cluster / group
#' (specified in \code{colnames.cluster}) and the selected variables is taken
#' and R squared values are computed based on the second halves of observations.
#' Finally, the R squared values are averaged over the \code{B} splits and over
#' the different data sets if multiple data sets are supplied.
#'
#' % Alternatively, second half-samples.
#'
#' For a continuous response, the adjusted R squared values is
#' calculated for a given cluster or group of variables. The Nagelkerke’s
#' R squared values is computed for a binary response using the function
#' \code{\link[fmsb]{NagelkerkeR2}}.
#'
#' If \code{colnames.cluster} is not supplied, the R squared value of the
#' full model is computed.
#'
#' @return The returned value is the R squared value.
#'
#' @seealso \code{\link[hierinf]{test_hierarchy}}.
#'
#' @examples
#' n <- 200
#' p <- 500
#' library(MASS)
#' set.seed(3)
#' x <- mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
#' colnames(x) <- paste0("Var", 1:p)
#' beta <- rep(0, p)
#' beta[c(5, 20, 46)] <- 1
#' y <- x %*% beta + rnorm(n)
#'
#' dendr <- cluster_var(x = x)
#' set.seed(47)
#' sign.clusters <- test_hierarchy(x = x, y = y, dendr = dendr,
#'                                 family = "gaussian")
#'
#' compute_r2(x = x, y = y, res.test.hierarchy = sign.clusters,
#'            family = "gaussian",
#'            colnames.cluster = c("Var1", "Var5", "Var8"))
#'
#' @references Renaux, C. et al. (2018), Hierarchical inference for genome-wide
#' association studies: a view on methodology with software. (arXiv:1805.02988)
#'
#' Nagelkerke, N. J. et al. (1991). A note on a general definition of the
#' coefficient of determination. Biometrika, 78:691–692.
#'
#' @name compute_r2
#' @export

compute_r2 <- function(x, y, res.test.hierarchy, clvar = NULL,
                       family = c("gaussian", "binomial"),
                       colnames.cluster = NULL) {

  family <- match.arg(family)

  if ("hierT" %in% class(res.test.hierarchy)) {
    res.multisplit <- res.test.hierarchy$res.multisplit
  } else {
    if ("hierM" %in% class(res.test.hierarchy)) {
      res.multisplit <- res.test.hierarchy
    } else {
      stop("The argument res.test.hierarchy is required to be the output of the function test_hierarchy, the function test_only_hierarchy, or the function multisplit.")
    }
  }

  # check input
  res <- check_input_testing(x = x, y = y, clvar = clvar, family = family,
                             # check result of the function multisplit
                             check_res_multisplit = TRUE,
                             res.multisplit = res.multisplit,
                             # arguments for the function multisplit
                             check_multisplit_arguments = FALSE,
                             B = NULL, proportion.select = NULL,
                             # arguments for the function
                             # test_hierarchy_given_multisplit
                             check_testing_arguments = FALSE,
                             dendr = NULL, block = NULL, alpha = NULL,
                             global.test = NULL, verbose = NULL)
  x <- res$x
  y <- res$y
  clvar <- res$clvar
  rm(list = c("res"))

  if (!is.null(attr(res.multisplit, "errorMsgs"))) {
    stop("There occured some errors in the previous function call of test_hierarchy or multisplit. Testing cannot be performed. See attribute 'errorMsgs' of the corresponding list element of the object which you specified in the argument res.test.hierarchy for more details.")
  }

  # Calculate unique.colnames.x
  len.x <- length(x) # this corresponds to the number of data sets
  colnames.x <- vector(mode = "character", length = 0)
  for (i in seq_len(len.x)) {
    colnames.x <- c(colnames.x, colnames(x[[i]]))
  }
  unique.colnames.x <- unique(x = colnames.x)

  # check colnames.cluster (it should only contains column names of x)
  if (!is.null(colnames.cluster)) {
    if (!all(colnames.cluster %in% unique.colnames.x)) {
      stop("Each variable which column names is specified in the argument colnames.cluster has to be contained in the data set or at least one data set for multiple data sets.")
    }
  }

  # Check that the selected variables in res.multisplit are contained in
  # unique.colnames.x
  colnames.multisplit <- vector(mode = "character", length = 0)
  for (i in seq_len(len.x)) {
    colnames.multisplit <- c(colnames.multisplit,
                             unique(as.vector(res.multisplit[[i]]$sel.coef)))
  }
  unique.coln.multisplit <- unique(x = colnames.multisplit)
  unique.coln.multisplit <- unique.coln.multisplit[!is.na(unique.coln.multisplit)]
  if (!all(unique.coln.multisplit %in% unique.colnames.x)) {
    stop("The selected variables in the output of the function call to test_hierarchy or multisplit does not match the column names of the argument x.")
  }

  # Defining the weights for aggregating the R^2 values
  weightR2 <- vapply(X = x, FUN = function(x) {nrow(x)}, FUN.VALUE = 1)
  weightR2 <- weightR2 / sum(weightR2)

  # compute the R2 value
  r2.data <- mapply(calculate_r2_one_data, x = x, y = y, clvar = clvar,
                    res.multisplit = res.multisplit,
                    MoreArgs = list(colnames.cluster = colnames.cluster,
                                    family = family))

  # compute the final R2
  return(sum(weightR2 * r2.data))
} # {compute_r2}

# Calculate R squared value for a given data set
#
# For a given data set, calculate R squared value for each of the \code{B}
# splits and takes the average
calculate_r2_one_data <- function(x, y, res.multisplit, clvar, family,
                                  colnames.cluster) {
  # prepare the variables for the call of comp_cluster_pval
  B <- nrow(res.multisplit$out.sample)

  # save all the rows of the matrix in a list
  out.sample <- split(res.multisplit$out.sample, seq(B))
  sel.coef <- split(res.multisplit$sel.coef, seq(B))

  # compute the p-value for each split and aggregate them
  r2.split <- mapply(FUN = calculate_r2_one_split, out.sample = out.sample,
                     sel.coef = sel.coef,
                     MoreArgs = list(x = x, y = y, clvar = clvar,
                                     colnames.cluster = colnames.cluster,
                                     family = family))

  return(mean(r2.split))
}

# Calculate R squared value for a given split and given data set
#
# For a given split and given data set, calculate the R squared value.
calculate_r2_one_split <- function(out.sample, sel.coef, x, y,
                                   clvar, family, colnames.cluster) {
  sel.coef <- sel.coef[!is.na(sel.coef)]

  if (is.null(colnames.cluster)) {
    # If colnames.cluster is equal to NULL, then calculate the R^2 of the entire
    # data set.
    common.colnames <- sel.coef
  } else {
    common.colnames <- intersect(colnames.cluster, sel.coef)
  }

  r2.one <-
    if (length(common.colnames) == 0) {
      0
    } else {
      return_r2(x = x[out.sample, common.colnames, drop = FALSE],
                y = y[out.sample], clvar = clvar[out.sample, ],
                family = family)
    }

  return(r2.one)

}

# Compute R squared value
#
# Calculate the R squared value for a linear or logistic regression model.
return_r2 <- function (x, y, clvar, family) {

  # generate design matrices
  design.mat <- cbind(clvar, x)
  # This results in a matrix although it might only have one column :-)
  # Note that if, say, clvar is equal to NULL, then this code works fine.
  # This means cbind(NULL, x) will result in x

  if (ncol(design.mat) == 0) {design.mat <- rep(1, length(y))}

  # compute r2
  r2 <-
    if (family == "binomial") {
      fmsb::NagelkerkeR2(MEL(design.mat, y, maxit = 100))$R2
    } else if (family == "gaussian") {
      stats::summary.lm(stats::lm(y ~ design.mat))$adj.r.squared
    }

  return(r2)
} # {return_r2}


