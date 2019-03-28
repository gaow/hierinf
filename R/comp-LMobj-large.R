# Compute model output for the large design matrix
#
# Compute model output for the large design matrix, i.e. based all the
# selected variables per split and per data set
#
# @return list of lm/glm objects for the large design matrix
compMOD_large <- function(x, y, clvar, res.multisplit, family) {

  MODobj_data <- mapply(compMOD_one_data, x = x, y = y, clvar = clvar,
                        res.multisplit = res.multisplit,
                        MoreArgs = list(family = family),
                        SIMPLIFY = FALSE)

  return(MODobj_data)
} # {compMOD_large}

# Compute output of lm/glm model for each data set
compMOD_one_data <- function(x, y, clvar, res.multisplit, family){

  # prepare the variables for the call of comp_cluster_pval
  B <- nrow(res.multisplit$out.sample)

  # save all the rows of the matrix in a list
  out.sample <- split(res.multisplit$out.sample, seq(B))
  sel.coef <- split(res.multisplit$sel.coef, seq(B))

  # compute the p-value for each split and aggregate them
  MODobj_split <- mapply(FUN = compMOD_one_split, out.sample = out.sample,
                         sel.coef = sel.coef,
                         MoreArgs = list(x = x, y = y, clvar = clvar,
                                         family = family),
                         SIMPLIFY = FALSE)

  return(MODobj_split)
} # {compMOD_one_data}

# Compute output of lm/glm model for each split
compMOD_one_split <- function(x, y, clvar, out.sample, sel.coef, family) {

  sel.coef <- sel.coef[!is.na(sel.coef)]

  MODobj <- compMOD_one(x = x[out.sample, sel.coef, drop = FALSE],
                        y = y[out.sample],
                        clvar = clvar[out.sample, ],
                        family = family)

  return(MODobj)
} # {compMOD_one_split}

# Compute output of lm/glm model

#' @importFrom stats lm
compMOD_one <- function (x, y, clvar, family) {

  # data.large <- cbind(clvar, x)

  # TODO use switch if there would be more possible families!
  MODout <-
    if (family == "binomial") {
      MEL(cbind(clvar, x), y, maxit = 100)
    } else if (family == "gaussian") {
      lm(y ~ cbind(clvar, x), model = FALSE, qr = FALSE)
    }

  return(MODout)
} # {compMOD_one}



