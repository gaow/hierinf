#' Multi-sample splitting
#'
#' The data is randomly split in two halves w.r.t. the observations and
#' variable selection using Lasso is performed on one half. Whereas the second
#' half and the selected variables are later used for testing by the function
#' \code{\link{test_only_hierarchy}}. This is repeated multiple times.
#'
#' @param x a matrix or list of matrices for multiple data sets. The matrix or
#' matrices have to be of type numeric and are required to have column names
#' / variable names. The rows and the columns represent the observations and
#' the variables, respectively.
#' @param y a vector, a matrix with one column, or list of the aforementioned
#' objects for multiple data sets. The vector, vectors, matrix, or matrices
#' have to be of type numeric. For \code{family = "binomial"}, the response
#' is required to be a binary vector taking values 0 and 1.
#' @param clvar a matrix or list of matrices of control variables.
#' @param B number of sample splits.
#' @param proportion.select proportion of variables to be selected by Lasso in
#' the multi-sample splitting step.
#' @param standardize a logical value indicating whether the variables should be
#' standardized.
#' @param family a character string naming a family of the error distribution;
#' either \code{"gaussian"} or \code{"binomial"}.
#' @param parallel type of parallel computation to be used. See the 'Details' section.
#' @param ncpus number of processes to be run in parallel.
#' @param cl an optional \strong{parallel} or \strong{snow} cluster used if
#' \code{parallel = "snow"}. If not supplied, a cluster on the local machine is created.
#' @param check.input a logical value indicating whether the function should
#' check the input. This argument is used to call
#' \code{\link{multisplit}} within
#' \code{\link{test_hierarchy}}.
#'
#' @details A given data with \code{nobs} is randomly split in two halves w.r.t.
#' the observations and \code{nobs * proportion.select} variables are selected
#' using Lasso (implemented in \code{\link{glmnet}}) on one half.
#' Control variables are not penalized if supplied
#' using the argument \code{clvar}. This is repeated \code{B} times for each
#' data set if multiple data sets are supplied. Those splits (i.e. second
#' halves of observations) and corresponding selected variables are used to
#' perform hierarchical testing by the function
#' \code{\link{test_only_hierarchy}}.
#'
#' The multi-sample split step can be run in parallel across the different
#' sample splits (\code{B} corresponds to number of sample splits) by
#' specifying the arguments \code{parallel} and \code{ncpus}.
#' There is an optional argument \code{cl} if \code{parallel = "snow"}.
#' There are three possibilities to set the argument \code{parallel}:
#' \code{parallel = "no"} for serial evaluation (default),
#' \code{parallel = "multicore"} for parallel evaluation
#' using forking, and \code{parallel = "snow"} for parallel evaluation
#' using a parallel socket cluster. It is recommended to select
#' \code{\link{RNGkind}("L'Ecuyer-CMRG")} and set a seed to ensure that
#' the parallel computing of the package \code{hierinf} is reproducible.
#' This way each processor gets a different substream of the pseudo random
#' number generator stream which makes the results reproducible if the arguments
#' (as \code{sort.parallel} and \code{ncpus}) remain unchanged. See the vignette
#' or the reference for more details.
#'
#' @return The returned value is an object of class \code{"hierM"}, consisting
#' of a list with number of elements corresponding to the number of data sets.
#' Each element (corresponding to a data set
#' % with \code{nobs} observations) contains a list with two matrices.
#' The first matrix
#' % of size \code{B x [nobs / 2]}
#' contains the indices of the second half of variables (which were not used
#' to select the variables). The second matrix
#' % of size \code{B x [nobs * proportion.select]}
#' contains the column names / variable names of the selected variables.
#'
#'
# A data frame with 2 components. A matrix of size \code{B x [nobs/2]}
# containing the second subsample of each split, and a matrix of size
# \code{B x [nobs/6]} containing the selected variables in each split.
#'
#' @seealso \code{\link{cluster_var}},
#' \code{\link{cluster_position}},
#' \code{\link{test_only_hierarchy}},
#' \code{\link{test_hierarchy}}, and
#' \code{\link{compute_r2}}.
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
#' set.seed(84)
#' res.multisplit <- multisplit(x = x, y = y, family = "gaussian")
#'
#' @references Renaux, C. et al. (2018), Hierarchical inference for genome-wide
#' association studies: a view on methodology with software. (arXiv:1805.02988)
#'
#' Meinshausen, N., Meier, L. and Buhlmann, P. (2009), P-values for
#' high-dimensional regression, Journal of the American Statistical Association
#' 104, 1671-1681.
#'
#' @name multisplit
#' @export

multisplit <- function(x, y, clvar = NULL, B = 50, proportion.select = 1/6,
                       standardize = FALSE, family = c("gaussian", "binomial"),
                       parallel = c("no", "multicore", "snow"), ncpus = 1L,
                       cl = NULL, check.input = TRUE) {

  family <- match.arg(family)
  parallel <- match.arg(parallel)
  do.parallel <- (parallel != "no" && ncpus > 1L)

  if (do.parallel && parallel == "multicore" && .Platform$OS.type == "windows") {
    stop("The argument parallel = 'multicore' is not available for windows. Use parallel = 'snow' for parallel computing or parallel = 'no' for serial execution of the code.")
  }

  ## check input
  if (check.input) {
    res <- check_input_testing(x = x, y = y, clvar = clvar, family = family,
                               # check result of the function multisplit
                               check_res_multisplit = FALSE,
                               res.multisplit = NULL,
                               # arguments for the function multisplit
                               check_multisplit_arguments = TRUE,
                               B = B, proportion.select = proportion.select,
                               standardize = standardize,
                               # arguments for the function
                               # test_hierarchy_given_multisplit
                               check_testing_arguments = FALSE,
                               dendr = NULL, block = NULL, alpha = NULL,
                               global.test = NULL, agg.method = NULL,
                               verbose = NULL)
    x <- res$x
    y <- res$y
    clvar <- res$clvar
    rm(list = c("res"))
  }

  ## create the splits
  len.x <- length(x)
  ret.split <- lapply(seq_len(len.x), create_splits_one_data, y = y, B = B)

  ## run the function multisplit_one_data for each of the data sets
  ret.sel <- lapply(seq_len(len.x), multisplit_one_data, x = x, y = y, clvar = clvar,
                    B = B, proportion.select = proportion.select,
                    standardize = standardize, family = family,
                    ret.split = ret.split, parallel = parallel,
                    ncpus = ncpus, cl = cl, do.parallel = do.parallel)

  RET <- do.call(cbind, ret.sel)
  resM.all <- structure(RET["resM", ], class = c("hierM", "list"))
  names(resM.all) <- NULL

  # do.call() returns NULL if there occurred no errors.
  attr(resM.all,"errorMsgs") <- do.call(c, RET["errorMsgs", ])
  attr(resM.all, "warningMsgs") <- do.call(c, RET["warningMsgs", ])

  if (!is.null(attr(resM.all, "errorMsgs"))) {
    warning("There occurred some errors while multi-sample splitting. See attribute 'errorMsgs' of the return object for more details.")
  }

  if (!is.null(attr(resM.all, "warningMsgs"))) {
    warning("There occurred some warnings while multi-sample splitting. See attribute 'warningMsgs' of the return object for more details.")
  }
  # return the values
  return(resM.all)
} # {multisplit}

# Create \code{B} splits for a given data set
#
# Returns indices for the \code{B} splits, i.e. in.sample and out.sample.
create_splits_one_data <- function(ind, y, B){

  # create the splits
  no.samples <- length(y[[ind]])

  in.sample <- lapply(seq_len(B), function(i, no.samples) {
    sort(sample(seq_len(no.samples), floor(no.samples/2)))
  }, no.samples = no.samples)

  out.sample <- lapply(in.sample, function(in.sample, no.samples) {
    setdiff(seq_len(no.samples), in.sample)
  }, no.samples = no.samples)

  return(list("in.sample" = in.sample, "out.sample" = out.sample))
} # {create_splits_one_data}

# Multi-sample split for a givne data set
#
# For each split, variables are selected using Lasso on the \code{in.sample}.
# A matrix of the selected variables (column names) and a matrix of
# indices of \code{out.sample} is returned. Both have \code{B} rows.
multisplit_one_data <- function(ind, x, y, clvar, B, proportion.select,
                                standardize, family, ret.split, parallel,
                                ncpus, cl, do.parallel){
  ## prepare variables
  in.out.sample <- ret.split[[ind]]

  y.one <- y[[ind]]
  no.samples <- length(y.one)

  xc <- cbind(clvar[[ind]], x[[ind]])
  no.xc <- ncol(xc)
  no.clvar <-
    if (is.null(clvar[[ind]])) {
      0
    } else {
      ncol(clvar[[ind]])
    }

  # don't penalize the control variables
  pen <- c(rep(0, no.clvar), rep(1, (no.xc - no.clvar)))

  # The concept of how to elegantly parallelize a function call (and save
  # all warning and error messages) is taken from the package boot
  # respectively lme4. Both are nearly identical in that respect.
  # See the source code of the package boot: R/bootfuns.q in the function
  # boot().
  # See the source code of the package lme4: R/bootMer.R in the function
  # bootMer().

  # Using a closure, the function below can access all the variables of the
  # environment in which it was created. This makes parallel computation
  # leaner or simpler, i.e. there are less arguments or we do not have to
  # export objects to the workers in the PSOCKcluster case
  one_split <- local({
    xc
    y.one
    ind
    no.clvar
    pen
    family
    proportion.select
    standardize
    no.samples
    in.out.sample
    function(i) {
      tryCatch_W_E(lasso_select(x = xc, y = y.one, no.clvar = no.clvar,
                                pen = pen, family = family,
                                proportion.select = proportion.select,
                                standardize = standardize,
                                no.samples = no.samples,
                                in.sample = in.out.sample$in.sample[[i]],
                                out.sample = in.out.sample$out.sample[[i]]),
                   ret.obj = list("out.sample" = rep(NA, no.samples - floor(no.samples/2)),
                                  "sel.coef" = rep(NA, no.samples * proportion.select)))
    }})

  ret <- if (do.parallel) {
    if (parallel == "multicore") {
      parallel::mclapply(seq_len(B), one_split, mc.cores = ncpus)
    } else if (parallel == "snow") {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        # export the namespace of hierINF in order for the use the functions
        # of the package hierINF on the workers
        parallel::clusterExport(cl, varlist = getNamespaceExports("hierINF"))
        if(RNGkind()[1L] == "L'Ecuyer-CMRG")
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, seq_len(B), one_split)
        parallel::stopCluster(cl)
        res
      } else parallel::parLapply(cl, seq_len(B), one_split)
    }
  } else lapply(seq_len(B), one_split)

  # Simplify the result
  ret <- do.call(cbind, ret)

  # prepare the return object
  rr <- do.call(cbind, ret["value", ])
  resM <- list("out.sample" = do.call(rbind, rr[1, ]),
               "sel.coef" = do.call(rbind, rr[2, ]))

  return(list("resM" = resM,
              "errorMsgs" = do.call(c, ret["error", ]),
              "warningMsgs" = do.call(c, ret["warning", ])))
} # {multisplit_one_data}

# Selects coefficents using Lasso
#
# For a given split, variables are selected using Lasso based on the
# observations / rows \code{in.sample}.
#
# @param pen a vector of 0 and 1 encoding the penalty. Control variables are
# unpenalized.
#
# @return a list with two elements, i.e. the indices of the second half of
# observations (\code{out.sample}) and the column / variable names of the
# selected variables (\code{sel.coef}).
lasso_select <- function(x, y, no.clvar, pen, family, proportion.select,
                         standardize, no.samples, in.sample, out.sample) {

  # perform lasso on the in.sample
  lasso.output <- glmnet::glmnet(x = x[in.sample, ], y = y[in.sample],
                                 family = family, standardize = standardize,
                                 intercept = TRUE, penalty.factor = pen)

  # save the first no.samples * proportion.select regression coefficients that
  # enter the lasso path
  lasso.coef.path <- lasso.output$beta@i + 1
  unique.lasso.coef <- unique(lasso.coef.path)

  # Note that no.samples * proportion.select may not be an integer but numeric
  # values are coerced to integer as by as.integer(), i.e. truncated towards zero.
  selected.lasso.coef <-
    if (no.clvar > 0) {
      setdiff(unique.lasso.coef, seq(1, no.clvar))[seq_len(no.samples * proportion.select)]
    } else {
      unique.lasso.coef[seq_len(no.samples * proportion.select)]
    }

  return(list("out.sample" = out.sample,
              "sel.coef" = rownames(lasso.output$beta)[selected.lasso.coef]))
} # {lasso_select}

