#' Hierarchical Testing
#'
#' Hierarchical Testing given the output of the function
#' \code{\link[hierinf]{multisplit}}.
#'
#' @param x a matrix or list of matrices for multiple data sets. The matrix or
#' matrices have to be of type numeric and are required to have column names
#' / variable names. The rows and the columns represent the observations and
#' the variables, respectively.
#' @param y a vector, a matrix with one column, or list of the aforementioned
#' objects for multiple data sets. The vector, vectors, matrix, or matrices
#' have to be of type numeric. For \code{family = "binomial"}, the response
#' is required to be a binary vector taking values 0 and 1.
#' @param dendr the output of one of the functions
#' \code{\link[hierinf]{cluster_var}} or \code{\link[hierinf]{cluster_position}}.
#' @param res.multisplit the output of the function
#' \code{\link[hierinf]{multisplit}}.
#' @param clvar a matrix or list of matrices of control variables.
#' @param family a character string naming a family of the error distribution;
#' either \code{"gaussian"} or \code{"binomial"}.
#' @param alpha the significant level at which the FWER is controlled.
#' @param global.test a logical value indicating whether the global test should
#' be performed.
#' @param verbose a logical value indicating whether the progress of the computation
#' should be printed in the console.
#' @param sort.parallel a logical indicating whether the values are sorted with respect to
#' the size of the block. This can reduce the run time for parallel computation.
#' @param parallel type of parallel computation to be used. See the 'Details' section.
#' @param ncpus number of processes to be run in parallel.
#' @param cl an optional \strong{parallel} or \strong{snow} cluster used if
#' \code{parallel = "snow"}. If not supplied, a cluster on the local machine is created.
#' @param check.input a logical value indicating whether the function should
#' check the input. This argument is used to call
#' \code{\link[hierinf]{test_only_hierarchy}} within
#' \code{\link[hierinf]{test_hierarchy}}.
#' @param unique.colnames.x a character vector containing the unique column
#' names of \code{x}. This argument is used to call
#' \code{\link[hierinf]{test_only_hierarchy}} within
#' \code{\link[hierinf]{test_hierarchy}}.
#'
#' @details The function \code{\link[hierinf]{test_only_hierarchy}} requires the output
#' of one of the functions \code{\link[hierinf]{cluster_var}} or
#' \code{\link[hierinf]{cluster_position}} as an input (argument \code{dendr}).
#' Furthermore it requires the output of the function
#' \code{\link[hierinf]{multisplit}} as an input (argument \code{res.multisplit}).
#' Hierarchical testing is performed by going top down through the hierarchical
#' tree. Testing only continues if at least one child of a given cluster is significant.
#'
#' If the argument \code{block} was supplied for the building
#' of the hierarchical tree (i.e. in the function call of either
#' \code{\link[hierinf]{cluster_var}} or
#' \code{\link[hierinf]{cluster_position}}), i.e. the second level of the
#' hierarchical tree was given, the hierarchical testing step can be run in
#' parallel across the different blocks by specifying the arguments
#' \code{parallel} and \code{ncpus}. There is an optional argument \code{cl} if
#' \code{parallel = "snow"}. There are three possibilities to set the
#' argument \code{parallel}: \code{parallel = "no"} for serial evaluation
#' (default), \code{parallel = "multicore"} for parallel evaluation
#' using forking, and \code{parallel = "snow"} for parallel evaluation
#' using a parallel socket cluster. It is recommended to select
#' \code{\link[base]{RNGkind}("L'Ecuyer-CMRG")} and set a seed to ensure that
#' the parallel computing of the package \code{hierinf} is reproducible.
#' This way each processor gets a different substream of the pseudo random
#' number generator stream which makes the results reproducible if the arguments
#' (as \code{sort.parallel} and \code{ncpus}) remain unchanged. See the vignette
#' or the reference for more details.
#'
#' @return The returned value is an object of class \code{"hierT"}, consisting of
#' two elements, the result of the multi-sample splitting step
#' \code{"res.multisplit"} and the result of the hierarchical testing
#' \code{"res.hierarchy"}.
#'
#' The result of the multi-sample splitting step is a list with number of
#' elements corresponding to the number of data sets. Each element
#' (corresponding to a data set) contains a list with two matrices. The first
#' matrix contains the indices of the second half of variables (which were
#' not used to select the variables). The second matrix contains the column
#' names / variable names of the selected variables.
#'
#' The result of the hierarchical testing is a data frame of significant
#' clusters with the following columns:
#' \item{block}{\code{NA} or the name of the block if the significant cluster
#' is a subcluster of the block or is the block itself.}
#' \item{p.value}{The p-value of the significant cluster.}
#' \item{significant.cluster}{The column names of the members of the significant
#' cluster.}
#'
#' There is a \code{print} method for this class; see
#' \code{\link[hierinf]{print.hierT}}.
#'
#' @seealso \code{\link[hierinf]{cluster_var}},
#' \code{\link[hierinf]{cluster_position}},
#' \code{\link[hierinf]{multisplit}},
#' \code{\link[hierinf]{test_hierarchy}}, and
#' \code{\link[hierinf]{compute_r2}}.
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
#' dendr1 <- cluster_var(x = x)
#' set.seed(76)
#' res.multisplit1 <- multisplit(x = x, y = y, family = "gaussian")
#' sign.clusters1 <- test_only_hierarchy(x = x, y = y, dendr = dendr1,
#'                                       res.multisplit = res.multisplit1,
#'                                       family = "gaussian")
#'
#' ## With block
#' # The column names of the data frame block are optional.
#' block <- data.frame("var.name" = paste0("Var", 1:p),
#'                     "block" = rep(c(1, 2), each = p/2),
#'                     stringsAsFactors = FALSE)
#' dendr2 <- cluster_var(x = x, block = block)
#' # The output res.multisplit1 can be used since the multi-sample
#' # step is the same with or without blocks.
#' sign.clusters2 <- test_only_hierarchy(x = x, y = y, dendr = dendr2,
#'                                       res.multisplit = res.multisplit1,
#'                                       family = "gaussian")
#'
#' # Access part of the object
#' sign.clusters2$res.hierarchy[, "block"]
#' sign.clusters2$res.hierarchy[, "p.value"]
#' # Column names or variable names of the significant cluster in the first row.
#' sign.clusters2$res.hierarchy[[1, "significant.cluster"]]
#'
#' @references Renaux, C. et al. (2018), Hierarchical inference for genome-wide
#' association studies: a view on methodology with software. (arXiv:1805.02988)
#'
#' @name test_only_hierarchy
#' @export

test_only_hierarchy <- function(x, y, dendr, res.multisplit, clvar = NULL,
                                family = c("gaussian", "binomial"),
                                alpha = 0.05, global.test = TRUE,
                                verbose = FALSE, sort.parallel = TRUE,
                                parallel = c("no", "multicore", "snow"),
                                ncpus = 1L, cl = NULL, check.input = TRUE,
                                unique.colnames.x = NULL) {

  block <- dendr$block
  dendr <- dendr$res.tree
  family <- match.arg(family)
  parallel <- match.arg(parallel)
  do.parallel <- (parallel != "no" && ncpus > 1L)

  if (do.parallel && parallel == "multicore" && .Platform$OS.type == "windows") {
    stop("The argument parallel = 'multicore' is not available for windows. Use parallel = 'snow' for parallel computing or parallel = 'no' for serial execution of the code.")
  }

  if (check.input) {
    res <- check_input_testing(x = x, y = y, clvar = clvar, family = family,
                               # check result of the function multisplit
                               check_res_multisplit = TRUE,
                               res.multisplit = res.multisplit,
                               # arguments for the function multisplit
                               check_multisplit_arguments = FALSE,
                               B = NULL, proportion.select = NULL,
                               standardize = NULL,
                               # arguments for the function
                               # test_hierarchy_given_multisplit
                               check_testing_arguments = TRUE,
                               dendr = dendr, block = block, alpha = alpha,
                               global.test = global.test, verbose = verbose)

    x <- res$x
    y <- res$y
    clvar <- res$clvar
    unique.colnames.x <- res$unique_colnames_x
    rm(list = c("res"))

    if (!is.null(attr(res.multisplit, "errorMsgs"))) {
      stop("There occurred some errors in the previous function call of multisplit. Testing cannot be performed. See attribute 'errorMsgs' of the object which you specified in the argument res.multisplit for more details.")
    }

    # Check that the selected variables in res.multisplit are contained in
    # unique.colnames.x
    colnames.multisplit <- lapply(X = res.multisplit,
                                  FUN = function(x) {
                                    unique(as.vector(x$sel.coef))
                                  })
    unique.coln.multisplit <- unique(unlist(colnames.multisplit))
    unique.coln.multisplit <- unique.coln.multisplit[!is.na(unique.coln.multisplit)]
    if (!all(unique.coln.multisplit %in% unique.colnames.x)) {
      stop("The selected variables in the output of the function call to test_hierarchy or multisplit does not match the column names of the argument x.")
    }
  }

  len.y <- length(y)
  if (verbose & len.y > 1) {
    message(paste("Jointly analyzing ", len.y, " phenotypes..."))
  }

  # Defining the weights for aggregating the p-values using Stouffer's method
  stouffer.weights <- vapply(X = x, FUN = function(x) {nrow(x)}, FUN.VALUE = 1)
  stouffer.weights <- sqrt(stouffer.weights / sum(stouffer.weights))

  # The variable minimal.pval is used for the hierarchical adjustment.
  # The p-value of a subcluster has to be as least as large as the p-value of
  # its parent.
  minimal.pval <- 0

  # This variable is used to stop testing if the global null hypothesis or all
  # the null hypotheses on the block level could not be rejected.
  continue.testing <- TRUE

  # This variable is used in order to store the warnings on the block level.
  warnings.to.return <- NULL

  ### testing the global null hypothesis ###
  if (global.test) {
    if (verbose) {
      message("Testing the global null hypothesis..")
    }
    # calculate the global p-value
    res.global <- tryCatch_W_E(comp_cluster_pval(x = x, y = y, clvar = clvar,
                                                 res.multisplit = res.multisplit,
                                                 colnames.cluster = unique.colnames.x,
                                                 family = family, len.y = len.y,
                                                 minimal.pval = minimal.pval,
                                                 stouffer.weights = stouffer.weights),
                               ret.obj = list("colnames.cluster" = NULL,
                                              "pval" = NULL))

    # If some warning occurred, then continue testing but report the warning
    # messages as an attribute of the return object.

    # If an error occurred during the computation of the global hypothesis,
    # then output all the error messages and stop running.
    if (!is.null(res.global$error)) {
      stop(paste("There occurred an errors while testing the global hypothesis.",
                 "All the error messages are printed below:",
                 paste(res.global$error,
                       collapse = "\n"),
                 "\n",
                 "All the warning messages are printed below:",
                 paste(res.global$warning,
                       collapse = "\n"),
                 sep = "\n"))
    }

    # Store warning messages. They are included as a attribute of the return
    # object.
    warnings.to.return <- res.global$warning

    # check if the global p-value is significant
    if (res.global$value$pval > alpha) {
      # the global p-value is larger than alpha
      if (verbose) {
        message("The global null hypothesis cannot be rejected.")
      }
      continue.testing <- FALSE
      signif.clusters <- list(list(value = list(name.block = NA,
                                                signif.clusters = list(
                                                  list(pval = NULL,
                                                       colnames.cluster = NULL))),
                                   error = NULL,
                                   warning = NULL))
    } else {
      # the global p-value is smaller than alpha => continue testing
      if (verbose) {
        message("The global null hypothesis was rejected.")
      }
      minimal.pval <- res.global$value$pval
    }
  }

  ### testing the blocks given by the argument block ###
  if (!is.null(block) & continue.testing) {
    if (verbose) {
      message("Testing the blocks...")
      # TODO find some better message: maybe subsets
      # testing .... number of blocks and their subsets
      # Testing the top clusters defined by the input block.
    }
    # test the blocks

    # The function split or divides the data x into blocks defined by f and stores
    # it in a list.
    colnames.per.block <- split(x = block[, 1], f = block[, 2])

    if (sort.parallel) {
      # Sort the blocks such that we test the large blocks first. This is
      # faster if we have less nodes / cpu's compared to the number of blocks.
      name.blocks <- names(sort(table(block[, 2]), decreasing = TRUE))
      colnames.per.block <- colnames.per.block[name.blocks]
    } else {
      name.blocks <- unique(block[, 2])
    }

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
    comp_per_blocks <- local({
      x
      y
      clvar
      res.multisplit
      family
      len.y
      minimal.pval
      stouffer.weights
      function(colnames.cluster) {
        tryCatch_W_E(comp_cluster_pval(x = x, y = y, clvar = clvar,
                                       res.multisplit = res.multisplit,
                                       colnames.cluster = colnames.cluster,
                                       family = family, len.y = len.y,
                                       minimal.pval = minimal.pval,
                                       stouffer.weights = stouffer.weights),
                     ret.obj = list("colnames.cluster" = NULL,
                                    "pval" = NULL))
      }})

    res.blocks <- if (do.parallel) {
      if (parallel == "multicore") {
        parallel::mclapply(colnames.per.block, comp_per_blocks, mc.cores = ncpus)
      } else if (parallel == "snow") {
        if (is.null(cl)) {
          cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
          # export the namespace of hierINF in order for the use the functions
          # of the package hierINF on the workers
          parallel::clusterExport(cl, varlist = getNamespaceExports("hierINF"))
          if(RNGkind()[1L] == "L'Ecuyer-CMRG")
            parallel::clusterSetRNGStream(cl)
          res <- parallel::parLapply(cl, colnames.per.block, comp_per_blocks)
          parallel::stopCluster(cl)
          res
        } else parallel::parLapply(cl, colnames.per.block, comp_per_blocks)
      }
    } else lapply(colnames.per.block, comp_per_blocks)

    res.blocks <- do.call(cbind, res.blocks)

    # If some warning occurred, then continue testing but report the warning
    # messages as an attribute of the return object.

    # If an error occurred during the computation per block, then output all the
    # error messages and stop running.
    if (!is.null(do.call(c, res.blocks["error", ]))) {
      stop(paste("There occurred one or more errors while testing each block.",
                 "All the error messages are printed below:",
                 paste(do.call(c, res.blocks["error", ]),
                       collapse = "\n"),
                 "\n",
                 "All the warning messages are printed below:",
                 paste(do.call(c, res.blocks["warning", ]),
                       collapse = "\n"),
                 sep = "\n"))
    }

    # Store warning messages. They are included as a attribute of the return
    # object.
    warnings.to.return <- c(warnings.to.return, do.call(c, res.blocks["warning", ]))

    # Check if any p-value of the blocks is significant.
    if (all(do.call(c, do.call(cbind, res.blocks["value", ])["pval", ]) > alpha)) {
      # All p-values of the blocks are larger than alpha
      if (verbose) {
        message("None of the null hypotheses for each block could be rejected.")
        message("Testing stops.")
      }
      if (global.test) {
        continue.testing <- FALSE
        signif.clusters <- list(list(value = list(name.block = NA,
                                                  signif.clusters =
                                                    list(res.global$value)),
                                     error = NULL,
                                     warning = NULL)) # See warnings.to.return
      } else {
        continue.testing <- FALSE
        signif.clusters <- list(list(value = list(name.block = NA,
                                                  signif.clusters = list(
                                                    list(pval = NULL,
                                                         colnames.cluster = NULL))),
                                     error = NULL,
                                     warning = NULL)) # See warnings.to.return
      }
    } else {
      # the p-value of the subset SNP_index is smaller than alpha => continue
      # testing
      if (verbose) {
        message("The null hypothesis of at least one block was rejected.")
        message("Testing the hierarchy of the corresponding significant blocks...")
        # TODO check the wording of the messages
      }
    }
  }

  ### prepare the input for the iterative testing for two special cases ###
  # Prepare the inputs for the function call of iterative_testing if the user
  # did not specify the argument block or the user did not specify the argument
  # block PLUS did set global.test to FALSE.
  if (is.null(block) & continue.testing) {
    # The function mapply cannot deal with arguments to vectorize over where
    # some arguments have strictly positive length and other arguments like
    # block have length 0. We use that list(NULL) has length 1 because it is
    # a list containing one element.
    name.blocks <- block

    # The second condition is to ensure that we test the top level of the tree
    # if the top level of the tree is not the same as the full data set /
    # global null. This makes it possible to use the package for each block,
    # say, chromosome separately as it is possible with the package hierGWAS.
    if (global.test) {
      if (length(setdiff(res.global$value$colnames.cluster, labels(dendr[[1]]))) == 0) {
        # top level of tree = global null
        test.top.level <- FALSE
      } else {
        # top level of tree != global null
        test.top.level <- TRUE
      }
    } else {
      # global.test = FALSE
      test.top.level <- TRUE
    }

    if (test.top.level) {
      # Top level of tree is tested
      #
      # There are two cases:
      # 1) global.test = FALSE: Test the top cluster of the tree in order to
      # initialize the iterative testing procedure. The result has to be
      # stored in a list with one element.
      # 2) global.test = TRUE: The top level of the tree is not the same as
      # the full data set / global null. Test first the top cluster of the
      # tree before continuing.
      res.blocks <- list(tryCatch_W_E(comp_cluster_pval(x = x, y = y,
                                                        clvar = clvar,
                                                        res.multisplit = res.multisplit,
                                                        colnames.cluster = labels(dendr[[1]]),
                                                        family = family,
                                                        len.y = len.y,
                                                        minimal.pval = minimal.pval,
                                                        stouffer.weights = stouffer.weights),
                                      ret.obj = list("colnames.cluster" = NULL,
                                                     "pval" = NULL)))

      res.blocks <- do.call(cbind, res.blocks)

      # If some warning occurred, then continue testing but report the warning
      # messages as an attribute of the return object.

      # If an error occurred during the computation of the top level of the tree (i.e. there are no blocks),
      # then output all the error messages and stop running.
      if (!is.null(do.call(c, res.blocks["error", ]))) {
        stop(paste("There occurred errors while testing the top level of the tree.",
                   "All the error messages are printed below:",
                   paste(do.call(c, res.blocks["error", ]),
                         collapse = "\n"),
                   "\n",
                   "All the warning messages are printed below:",
                   paste(do.call(c, res.blocks["warning", ]),
                         collapse = "\n"),
                   sep = "\n"))
      }

      # Store warning messages. They are included as a attribute of the return
      # object.
      warnings.to.return <- c(warnings.to.return, do.call(c, res.blocks["warning", ]))

      # If global.test = TRUE and the top cluster of the tree is not
      # significant, then return the res.global. (If global.test = FALSE, then
      # the function iterative_testing takes care of that.)
      # We could ommit the all() below because it's just one value but that
      # does not hurt.
      if (global.test & all(do.call(c, do.call(cbind, res.blocks["value", ])["pval", ]) > alpha)) {
        # All p-values of the blocks are larger than alpha
        if (verbose) {
          message("The null hypotheses of the top level of the tree could not be rejected.")
          message("Testing stops.")
        }
        continue.testing <- FALSE
        signif.clusters <- list(list(value = list(name.block = NA,
                                                  signif.clusters =
                                                    list(res.global$value)),
                                     error = NULL,
                                     warning = NULL)) # See warnings.to.return
      }
    } else {
      # top level of tree does not have to be tested
        res.blocks <- list(res.global)

        res.blocks <- do.call(cbind, res.blocks)
    }
  }

  ### testing the hierarchy defined by the tree (for all significant blocks) ###
  if (continue.testing) {
    # Sort the list of dendrograms.
    # It is needed if, say, sort.parallel is set to FALSE for the building of
    # the hierarchical tree but is set to TRUE for the function call of
    # test_hierarchy.
    if (!is.null(block)) {
      dendr <- dendr[name.blocks]
    }

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
    cluster_the_blocks <- local({
      x
      y
      clvar
      res.multisplit
      family
      len.y
      alpha
      verbose
      dendr
      name.blocks
      res.blocks
      stouffer.weights
      function(i) {
        tryCatch_W_E(iterative_testing(x = x, y = y, clvar = clvar,
                                       res.multisplit = res.multisplit,
                                       dendr = dendr[[i]],
                                       name.block = name.blocks[i],
                                       res.block = res.blocks[["value", i]],
                                       family = family,
                                       len.y = len.y, alpha = alpha,
                                       verbose = verbose,
                                       stouffer.weights = stouffer.weights),
                     ret.obj = list(name.block = name.blocks[i],
                                    signif.clusters = list(list(colnames.cluster = NULL,
                                                                pval = NULL))))
      }})


    # The sorting is done above during the testing of the block level.
    ind <- seq_len(dim(res.blocks)[2])

    signif.clusters <- if (do.parallel) {
      if (parallel == "multicore") {
        parallel::mclapply(ind, cluster_the_blocks, mc.cores = ncpus)
      } else if (parallel == "snow") {
        if (is.null(cl)) {
          cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
          # export the namespace of hierINF in order for the use the functions
          # of the package hierINF on the workers
          parallel::clusterExport(cl, varlist = getNamespaceExports("hierINF"))
          if(RNGkind()[1L] == "L'Ecuyer-CMRG")
            parallel::clusterSetRNGStream(cl)
          res <- parallel::parLapply(cl, ind, cluster_the_blocks)
          parallel::stopCluster(cl)
          res
        } else parallel::parLapply(cl, ind, cluster_the_blocks)
      }
    } else lapply(ind, cluster_the_blocks)
  }

  signif.clusters <- do.call(cbind, signif.clusters)
  sig.cl.compact <- lapply(X = signif.clusters["value", ], FUN = prepare_output)
  sig.cl.compact <- do.call(rbind, sig.cl.compact)
  colnames(sig.cl.compact) <- c("block", "p.value", "significant.cluster")
  rownames(sig.cl.compact) <- NULL

  resT <- sig.cl.compact
  # do.call() returns NULL if there occurred no errors.
  attr(resT,"errorMsgs") <- do.call(c, signif.clusters["error", ])
  # Add warning messages from block levels
  attr(resT, "warningMsgs") <- c(warnings.to.return,
                                 do.call(c, signif.clusters["warning", ]))

  if (!is.null(attr(resT, "errorMsgs"))) {
    warning("There occurred some errors while testing the hierarchy. See attribute 'errorMsgs' of the corresponding list element of the return object for more details.")
  }

  if (!is.null(attr(resT, "warningMsgs"))) {
    warning("There occurred some warnings while testing the hierarchy. See attribute 'warningMsgs' of the corresponding list element of the return object for more details.")
  }

  retH <- list("res.multisplit" = res.multisplit, "res.hierarchy" = resT)
  retH <- structure(retH, class = c("hierT", "list"))

  return(retH)
} # {test_hierarchy}

# Prepare the output
#
# This function changes the format of the output.
prepare_output <- function(signif.clusters) {
  name.block <-
    if (is.null(signif.clusters$name.block)) {
      NA
    } else {
      signif.clusters$name.block
    }
  len.sig.cl <- length(signif.clusters$signif.clusters)
  res.out <- cbind(name.block, data.frame(matrix(NA, nrow = len.sig.cl, ncol = 2)),
                   stringsAsFactors = FALSE)
  for (i in seq_len(len.sig.cl)) {
    res.out[i, 2:3] <-
      if (!is.null(signif.clusters$signif.clusters[[i]]$colnames.cluster)) {
        # We use a list instead of a vector in order for the p-value not to be
        # converted to a character.
        list(signif.clusters$signif.clusters[[i]]$pval,
             list(signif.clusters$signif.clusters[[i]]$colnames.cluster))
      } else {
        list(NA, list(NA)) # list(NA, NA) would work as well.
      }
  }
  return(res.out)
} # {prepare_output}

