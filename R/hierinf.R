#' hierinf: Hierarchical Inference
#'
#' The hierinf package provides the functions to perform hierarchical inference.
#' The main workflow consists of two function calls.
#'
#' @section Functions:
#' The building of the hierarchical tree can be achieved by either of the functions
#' \code{\link[hierinf]{cluster_var}} or \code{\link[hierinf]{cluster_position}}.
#' The function \code{\link[hierinf]{test_hierarchy}} performs the hierarchical
#' testing by going top down through the hierarchical tree and obviously requires the
#' hierarchical tree as an input.
#'
#' It is possible to calculate the R squared value of a given cluster using
#' \code{\link[hierinf]{compute_r2}}.
#'
#' The hierarchical testing consists of two steps which can be evaluated
#' separately if desired. Instead of calling \code{\link[hierinf]{test_hierarchy}},
#' the multi-sample splitting step is performed by \code{\link[hierinf]{multisplit}}
#' and its output is used by the function \code{\link[hierinf]{test_only_hierarchy}}
#' to test for significant clusters by going top down through the hierarchical tree.
#'
#' @docType package
#' @name hierinf
NULL
