# Compute the adjusted p-value of a given cluster
#
# Compute the adjusted p-value of a given cluster (specified by the argument
# \code{colnames.cluster}). This means that there is one adjusted p-value
# based on all data sets if multiple data sets are supplied. The p-values
# per data set are combined using Stouffer's method.
#
# @return adjusted p-value.
comp_cluster_pval <- function(x, y, clvar, res.multisplit, colnames.cluster,
                              family, len.y, minimal.pval, stouffer.weights) {

  # compute a p-value for each of the phenotypes or phenotypes & corresponding
  # (distinct) genotypes and (distinct) control covariates
  pvals_data <- mapply(comp_one_data, x = x, y = y, clvar = clvar,
                       res.multisplit = res.multisplit,
                       MoreArgs = list(colnames.cluster = colnames.cluster,
                                       family = family))

  # Stouffer's method: combine the p-values
  pval <- stats::pnorm(sum(stouffer.weights * stats::qnorm(pvals_data)))

  # hierarchical adjustment of the p-value (below Equation 4 on page 333 of
  # Mandozzi and Buehlmann (2016))
  return(list(colnames.cluster = colnames.cluster, pval = max(pval, minimal.pval)))
} # {comp_cluster_pval}

# Compute the adjusted p-value for a given cluster and given data set
#
# Compute the adjusted p-value for a given cluster (specified by the
# argument \code{colnames.cluster}) and given data set.
comp_one_data <- function(x, y, clvar, res.multisplit, colnames.cluster, family){

  # prepare the variables for the call of comp_cluster_pval
  B <- nrow(res.multisplit$out.sample)

  # save all the rows of the matrix in a list
  out.sample <- split(res.multisplit$out.sample, seq(B))
  sel.coef <- split(res.multisplit$sel.coef, seq(B))

  # compute the p-value for each split and aggregate them
  pvals.split <- mapply(FUN = comp_one_split, out.sample = out.sample,
                        sel.coef = sel.coef,
                        MoreArgs = list(x = x, y = y, clvar = clvar,
                                        colnames.cluster = colnames.cluster,
                                        family = family))
  if ((no_NA <- sum(is.na(pvals.split))) == B) {
    stop("The p-value of a cluster could not be calculated for all of the ", B, " splits. You might have colinear variables in one of your data sets and the algorithm might try to test a cluster containing only colinear variables. There could be more than one colinear variable in that cluster but not all of them.")
  }

  if (no_NA > 0) {
    warning(paste0("The p-value of a cluster could not be calculated in ", no_NA, " out of ", B, " splits for one of the data sets. This might be a problem because of colinear variables. The algorithm might try to test a cluster containing only colinear variables (but not all of them)."))
  }

  # Aggregation of p-values per split
  # Equation 4 on page 333 in Mandozzi and Buehlmann (2016)
  return(adj_pval(pvals.split, B))
} # {comp_one_data}

# Compute the adjusted p-value for a given cluster and given split of a data
# set
#
# Compute the adjusted p-value for a given cluster (specified by the
# argument \code{colnames.cluster}) and given split of a data set.
comp_one_split <- function(x, y, clvar, out.sample, sel.coef, colnames.cluster,
                           family) {
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
                             family = family)
      # Equation 3 on page 333 in Mandozzi and Buehlmann (2016)
      min(pval_unadj * length(sel.coef) / length(common.colnames), 1)
    }
  # return adjusted p-value
  return(pval)
} # {comp_one_split}

# Perform LRT
#
# Perform LRT (or F test) and return the resulting p-value.
test_var <- function (x, y, clvar, colnames.cluster, family) {

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
      stats::anova(MEL(data.small, y, maxit = 100),
                   MEL(data.large, y, maxit = 100),
                   test = "Chisq")$"Pr(>Chi)"[2]
    } else if (family == "gaussian") {
      # partial F test
      stats::anova(stats::lm(y ~ data.small),
                   stats::lm(y ~ data.large),
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

