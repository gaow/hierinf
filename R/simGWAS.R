#' Simulated GWAS data set
#'
#' The data set \code{simGWAS} was simulated using \code{PLINK} where the SNPs
#' were binned into different allele frequency ranges. There are 250 controls
#' and 250 cases, i.e. a binary response and 500 subjects.
#' The variables \code{age} and \code{sex}
#' are two additional control variables. The variables \code{SNP.1} till
#' \code{SNP.990} were simulated to have no association with the response and
#' the variables \code{SNP.991} till \code{SNP.1000} have a population odds
#' ratio of 2.
#'
#' @usage data(simGWAS)
#'
#' @format A list with three elements:
#' \describe{
#'    \item{\code{x}}{a matrix with 500 rows and 1000 columns where the rows
#'    and columns correspond to the subjects and variables, respectively.
#'    The variables are named \code{SNP.1}, ..., \code{SNP.1000}.}
#'    \item{\code{y}}{binary response vector with 500 elements where the
#'    elements correspond to the subjects.}
#'    \item{\code{clvar}}{a matrix with 500 rows and 2 columns where the rows
#'    and columns correspond to the subjects and variables, respectively.
#'    The age of the subject is stored in the variable \code{age}. The
#'    variable \code{sex} takes the value 0 for men and 1 for women.}
#' }
#'
#' @source
#' Buzdugan L (2018). hierGWAS: Asessing statistical significance in predictive
#' GWA studies. R package version 1.10.0.
#'
#' @examples
#' data(simGWAS)
#' sim.geno <- simGWAS$x
#' sim.pheno <- simGWAS$y
#' sim.clvar <- simGWAS$clvar
#'
#' dendr <- cluster_var(x = sim.geno)
#' result <- test_hierarchy(x = sim.geno, y = sim.pheno,
#'                          dendr = dendr, clvar = sim.clvar,
#'                          family = "binomial", seed = 1234)
"simGWAS"
