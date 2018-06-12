# Iterative testing
#
# The function goes top down through a given hierarchical tree and returns
# the significant clusters or single variables. A cluster is returned as
# a significant cluster if the cluster is significant and both children
# are not significant (if the cluster isn't a leaf). The list
# \code{cluster.queue} keeps track of the cluster to be tested.
iterative_testing <- function(x, y, clvar, res.multisplit, dendr, name.block,
                              res.block, family, len.y, alpha, verbose,
                              stouffer.weights) {

  ### check if the current cluster is significant ###
  if (res.block$pval > alpha){
    return(list(name.block = name.block,
                signif.clusters = list(list(colnames.cluster = NULL,
                                            pval = NULL))))
  }

    # Initialize the cluster.queue
  # This is a list with one element which contains a list with elements.
  cluster.queue <- list(list(dendr = dendr,
                             colnames.cluster = res.block$colnames.cluster,
                             pval = res.block$pval))
  # Note that res.block$pval is the same as minimal.pval for the children of
  # that cluster.

  # TODO maybe remove this line
  stopifnot(length(setdiff(labels(dendr), res.block$colnames.cluster)) == 0)

  # the significant clusters are stored in this following list
  signif.clusters <- vector("list", 0)

  # the whole functions is a big loop over the list cluster.queue
  while (length(cluster.queue) > 0) {
    ### get the current cluster ###

    # extract some values of the first entry of the list
    clust <- cluster.queue[[1]]$dendr
    size.cluster <- length(cluster.queue[[1]]$colnames.cluster)

    if (size.cluster == 1) {
      ### check if the current cluster is a leaf ###

      # put the significant leaf in the list of results
      signif.clusters[[length(signif.clusters) + 1]] <-
        list(pval = cluster.queue[[1]]$pval,
             colnames.cluster = cluster.queue[[1]]$colnames.cluster)

      # remove first element of the list cluster.queue
      cluster.queue <- cluster.queue[-1]
    } else {
      ### The current cluster is not a leaf => continue testing ###
      clust.height <- attr(clust, "height")

      if (verbose & !is.null(name.block)) {
        message(paste0("For block ", name.block, " testing the children of a cluster at height ", clust.height, "."))
      } else if (verbose & is.null(name.block)) {
        message(paste0("Testing the children of a cluster at height ", clust.height, "."))
      }
      # TODO check the wording of the messages

      # get the subclusters of the current cluster
      subclust <- cut(clust, h = clust.height)$lower
      fun_labels <- function(x, subclust) {
        labels(subclust[[x]])
      }
      colnames.subcluster <- lapply(seq_along(subclust), fun_labels,
                                    subclust = subclust)
      ### check the children of this cluster ###
      res.children <- mapply(FUN = comp_cluster_pval,
                             colnames.cluster = colnames.subcluster,
                             MoreArgs = list(x = x, y = y,
                                             clvar = clvar,
                                             res.multisplit = res.multisplit,
                                             family = family, len.y = len.y,
                                             minimal.pval = cluster.queue[[1]]$pval,
                                             stouffer.weights = stouffer.weights),
                             SIMPLIFY = FALSE)

      # check which clusters are significant
      update.signif.clusters <- mapply(FUN = check_significant,
                                       res.child = res.children,
                                       subcluster = subclust,
                                       MoreArgs = list(alpha = alpha),
                                       SIMPLIFY = FALSE)
      ind.NULL <- vapply(update.signif.clusters,
                         FUN = function(x) {ifelse(is.null(x), TRUE, FALSE)},
                         FUN.VALUE = TRUE)

      if (sum(!ind.NULL) > 0){
        # Add the significant children to the cluster.queue
        cluster.queue <- c(cluster.queue, update.signif.clusters[!ind.NULL])
      } else {
        # None of the children is significant
        # put the significant leaf in the list of results
        signif.clusters[[length(signif.clusters) + 1]] <-
          list(pval = cluster.queue[[1]]$pval,
               colnames.cluster = cluster.queue[[1]]$colnames.cluster)
      }

       # remove current cluster from cluster.queue
       cluster.queue <- cluster.queue[-1]
    }
  }

  if (verbose & !is.null(name.block)) {
    message(paste0("There have been found ", length(signif.clusters), " significant clusters for block ", name.block, ". Stop testing block ", name.block, "."))
  } else if (verbose & is.null(name.block)) {
    message(paste0("There have been found ", length(signif.clusters), " significant clusters. Stop testing."))
  }
  # TODO check the wording of the messages

  # return the list of significant clusters
  return(list(name.block = name.block,
              signif.clusters = signif.clusters))
} # {iterative_testing}

# Check if a cluster (it's a child of some cluster) is significant
#
# Check if a cluster (it's a child of some cluster) is significant. It
# is significant, then it has to be added to the list \code{cluster.queue}
# inside the function \code{iterative_testing}.
check_significant <- function(res.child, alpha, subcluster) {
  ret <-
    if(res.child$pval > alpha) {
      NULL
    } else {
      list(dendr = subcluster,
           colnames.cluster = res.child$colnames.cluster,
           pval = res.child$pval)
    }
  return(ret)
} # {check_significant}

