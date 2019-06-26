# Compute the adjusted p-value of a given cluster
#
# Compute the adjusted p-value of a given cluster (specified by the argument
# \code{colnames.cluster}). This means that there is one adjusted p-value
# based on all data sets if multiple data sets are supplied. The p-values
# per data set are combined using Stouffer's method.
#
# @return adjusted p-value.
comp_cluster_pval <- function(x, y, clvar, res.multisplit, colnames.cluster,
                              family, len.y, minimal.pval, agg.method, mod.large,
                              stouffer.weights) {

  # compute a p-value for each of the phenotypes or phenotypes & corresponding
  # (distinct) genotypes and (distinct) control covariates
  pvals_data <- mapply(comp_one_data, x = x, y = y, clvar = clvar,
                       res.multisplit = res.multisplit, mod.large = mod.large,
                       MoreArgs = list(colnames.cluster = colnames.cluster,
                                       family = family))

  pval <-
    if (length(pvals_data) == 1) {
      # No aggregation method is applied because the user only specified one
      # data set.
      pvals_data

    } else if (agg.method == "Tippett") {
      # Tippett's method: combine the p-values
      max(1 - (1 - min(pvals_data))^(length(x)), .Machine$double.neg.eps)
      # We use max(., .Machine$double.neg.eps) because all smaller values
      # are set to zero, i.e. identical(1, 1 - 1e-17) => TRUE because of
      # rounding in floating point arithmetic.

      # # Alternative:
      # # Minimum p-value, Bonferroni corrected, i.e. m * min(p_i)
      # min(c(1, length(x) * min(pvals_data)))

    } else if (agg.method == "Stouffer") {
      # Stouffer's method: combine the p-values
      stats::pnorm(sum(stouffer.weights * stats::qnorm(pvals_data)))

    }
    # else if (agg.method == "max") {
    #   # Largest p-value
    #   max(pvals_data)^(length(x))
    # }


  # hierarchical adjustment of the p-value (below Equation 4 on page 333 of
  # Mandozzi and Buehlmann (2016))
  return(list(colnames.cluster = colnames.cluster, pval = max(pval, minimal.pval)))
} # {comp_cluster_pval}

# Compute the adjusted p-value for a given cluster and given data set
#
# Compute the adjusted p-value for a given cluster (specified by the
# argument \code{colnames.cluster}) and given data set.
comp_one_data <- function(x, y, clvar, res.multisplit, colnames.cluster,
                          family, mod.large){

  # prepare the variables for the call of comp_cluster_pval
  B <- nrow(res.multisplit$out.sample)

  # save all the rows of the matrix in a list
  out.sample <- split(res.multisplit$out.sample, seq(B))
  sel.coef <- split(res.multisplit$sel.coef, seq(B))

  # compute the p-value for each split and aggregate them
  pvals.split <- mapply(FUN = comp_one_split, out.sample = out.sample,
                        sel.coef = sel.coef, mod.large = mod.large,
                        MoreArgs = list(x = x, y = y, clvar = clvar,
                                        colnames.cluster = colnames.cluster,
                                        family = family))

  if ((no_NA <- sum(is.na(pvals.split))) > 0) {
    warning(paste0("The p-value of a cluster could not be calculated in ", no_NA, " out of ", B, " splits for one of the data sets. The corresponding p-values are set to 1. This problem can occure due to colinear variables which can be linear combinations of variables. The algorithm might try to test a cluster containing (only) colinear variables but not all of them."))
    pvals.split[is.na(pvals.split)] <- 1
  }

  # Aggregation of p-values over the B splits
  # Equation 4 on page 333 in Mandozzi and Buehlmann (2016)
  return(adj_pval(pvals.split, B))
} # {comp_one_data}

# Compute the adjusted p-value for a given cluster and given split of a data
# set
#
# Compute the adjusted p-value for a given cluster (specified by the
# argument \code{colnames.cluster}) and given split of a data set.
comp_one_split <- function(x, y, clvar, out.sample, sel.coef, colnames.cluster,
                           family, mod.large) {
  sel.coef <- sel.coef[!is.na(sel.coef)]
  common.colnames <- intersect(colnames.cluster, sel.coef)

  # maybe change this !
  pval <-
    if (length(common.colnames) == 0) {
      1 # The p-value does not have to be calculated.
    } else {
      # drop = FALSE because we need a matrix although only one column might be
      # selected.
      pval_unadj <- test_var(x = x[out.sample, sel.coef, drop = FALSE],
                             y = y[out.sample],
                             clvar = clvar[out.sample, ],
                             colnames.cluster = colnames.cluster,
                             family = family,
                             mod.large = mod.large)
      # Equation 3 on page 333 in Mandozzi and Buehlmann (2016)
      min(pval_unadj * length(sel.coef) / length(common.colnames), 1)
    }
  # return adjusted p-value
  return(pval)
} # {comp_one_split}

# Perform LRT
#
# Perform LRT (or F test) and return the resulting p-value.

#' @importFrom stats lm anova
test_var <- function (x, y, clvar, colnames.cluster, family, mod.large) {

  ### generate design matrices ###
  setdiff.cluster <- setdiff(colnames(x), colnames.cluster)

  data.large <- cbind(clvar, x)
  data.small <- cbind(clvar, x[, setdiff.cluster]) # This results in a matrix although it might only have one column :-)
  # Note that if, say, clvar is equal to NULL, then this code works fine.
  # This means cbind(NULL, x) will result in x

  ### compare the models ###
  if (ncol(data.small) == 0) {data.small <- rep(1, length(y))}
  # TODO use switch if there would be more possible families!
  pval <-
    if (family == "binomial") {
      # likelihood ratio test

      # stats::anova(MEL(data.small, y, maxit = 100),
      #              # MEL(data.large, y, maxit = 100),
      #              mod.large,
      #              test = "Chisq")$"Pr(>Chi)"[2]

      own_anova.glmlist(list(MEL(data.small, y, maxit = 100),
                             # MEL(data.large, y, maxit = 100),
                             mod.large),
                        test = "Chisq")$"Pr(>Chi)"[2]
    } else if (family == "gaussian") {
      # partial F test
      anova(lm(y ~ data.small, model = FALSE, qr = FALSE),
            # stats::lm(y ~ data.large),
            mod.large,
            test = "F")$P[2]
    }

  return(pval)
} # {test_var}

# Adjust and aggregate the p-values (per split)
#
# Adjust and aggregate the \code{B} p-values (per split) for a given cluster
# and given data set.
adj_pval <- function(pvals, B) {
  # define the sequence of gamma values
  gamma.min <- 0.05
  gamma.step <- 0.01
  gamma.seq <- seq(gamma.min, 1, gamma.step)

  # compute the empirical quantile vector
  gamma.step <- vapply(X = gamma.seq,
                       FUN = function(g, pvals) {
                         min(1, stats::quantile(pvals / g, g, na.rm = TRUE))
                       },
                       FUN.VALUE = numeric(1),
                       pvals = pvals)

  # compute the adjusted p value
  # Equation 4 on page 333 in Mandozzi and Buehlmann (2016)
  return(min(1, (1 - log(gamma.min)) * min(gamma.step)))
} # {adj_pval}

