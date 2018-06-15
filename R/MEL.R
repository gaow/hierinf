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

MEL <- function(x, y, maxit, delta = 0.01, epsilon = 1e-6) {
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
                                maxit = maxit)))

  return(res)
}



