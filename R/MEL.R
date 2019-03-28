# Algorithm published in:
# Rousseeuw, Peter J., and Andreas Christmann. Robustness against separation
# and outliers in logistic regression. Computational Statistics & Data Analysis
# 43.3 (2003): 315-332.
#
# Instead of having observations with values 0 and 1, we create
# pseudo-observations of the form, say, c(0.04, 0.96) whereas with group data
# the first number corresponds to the "number" of successes and the second
# number corresponds to the "number" of failures. The advantage is that
# we do not have the problem of separation or quasi-separation of logistic
# regression since the maximum estimated likelihood (MEL) estimator always
# exists compared to the classical maximum likelihood.
#
# They recommend to use delta equal to 0.01.


# We import the functions glm and glm.control from the stats package becuase
# those functions are called many times for family = "binomial".
# See book http://r-pkgs.had.co.nz/namespace.html
# Wickham, Hadley. R packages: organize, test, document, and share your code.
# O'Reilly Media Inc., 2015.
# "If you are using functions repeatedly, you can avoid :: by importing the
# function with @importFrom pkg fun. This also has a small performance benefit,
# because :: adds approximately 5 microsecond to function evaluation time."

#' @importFrom stats glm glm.control

MEL <- function(x, y, maxit, delta = 0.01, epsilon = 1e-6, model = FALSE) {
  # mean of y but we bound it awy from 0 and 1. See Equation (10) on page 8.
  pi.hat <- max(delta, min(1 - delta, mean(y)))
  # The two parameters delta.0 and delta.1 are constrained such that the average
  # of the pseudo-observation is equal to pi.hat. See Equation (9) on page 8.
  delta.0 <- (pi.hat * delta) / (1 + delta)
  delta.1 <- (1 + pi.hat * delta) / (1 + delta)
  # Pseudo-observations. See Equation (3) on page 4.
  y.tilde <- delta.0 * (1 - y) + delta.1 * y
  pseudo.y <- cbind(y.tilde, 1 - y.tilde)

  # Suppress warning that "non-integer counts in a binomial glm!" because the
  # function glm expects in the first column to be the number of successes and
  # the second column to be the number of failures
  suppressWarnings(res <- glm(pseudo.y ~ x, family = "binomial",
                              control = glm.control(
                                epsilon = epsilon,
                                maxit = maxit),
                              model = model, y = FALSE))

  # Two parts of the output object to decrease the size of the object.
  res$R <- NULL
  res$qr$qr <- NULL
  # res$qr$qr <- matrix(NA, ncol = 2, nrow = 2)

  return(res)
}



# MODIFIED version of stats:::anova.glmlist in order to be able to store
# the output of our function MEL / the output of glm more compact, i.e.
# this function allows to drop the element "qr" of the output of glm
# because the function own_anova.glmlist does not call summary.

#' @importFrom stats formula stat.anova
own_anova.glmlist <- function (object, dispersion = NULL, test = NULL) {
  ### we could drop this part ###
  responses <- as.character(lapply(object, function(x) {
    deparse(formula(x)[[2L]])
  }))
  sameresp <- responses == responses[1L]
  if (!all(sameresp)) {
    object <- object[sameresp]
    warning(gettextf("models with response %s removed because response differs from model 1",
                     sQuote(deparse(responses[!sameresp]))), domain = NA)
  }
  ns <- sapply(object, function(x) length(x$residuals))
  if (any(ns != ns[1L]))
    stop("models were not all fitted to the same size of dataset")
  ### End: we could drop this part ###

  # [...]

  nmodels <- length(object)
  # if (nmodels == 1)
  #   return(anova.glm(object[[1L]], dispersion = dispersion,
  #                    test = test))
  resdf <- as.numeric(lapply(object, function(x) x$df.residual))
  resdev <- as.numeric(lapply(object, function(x) x$deviance))
  # [...]
  table <- data.frame(resdf, resdev, c(NA, -diff(resdf)),
                      c(NA, -diff(resdev)))
  variables <- lapply(object, function(x) paste(deparse(formula(x)),
                                                collapse = "\n"))
  dimnames(table) <- list(1L:nmodels, c("Resid. Df", "Resid. Dev",
                                        "Df", "Deviance"))
  # [...]
  title <- "Analysis of Deviance Table\n"
  topnote <- paste0("Model ", format(1L:nmodels), ": ", variables,
                    collapse = "\n")
  if (!is.null(test)) {
    bigmodel <- object[[order(resdf)[1L]]]
    dispersion <- 1 # MODIFIED
    df.dispersion <- Inf # MODIFIED
    table <- stat.anova(table = table, test = test, scale = dispersion,
                        df.scale = df.dispersion,
                        n = length(bigmodel$residuals))
  }
  structure(table, heading = c(title, topnote),
            class = c("anova", "data.frame"))
}



