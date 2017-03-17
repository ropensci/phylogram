#' Reorder tree branches in ladderized format.
#'
#' \code{ladder} takes an object of class \code{dendrogram} and ladderizes
#'   the branches for easier visual interpretation.
#'
#' @param x an object of class \code{"dendrogram"}.
#' @param decreasing logical indicating whether the tree should be
#'   ladderized upwards (FALSE; default) or downwards (TRUE).
#' @return an object of class \code{dendrogram}.
#' @details TBA.
#' @author Shaun Wilkinson
#' @seealso The \code{\link[ape]{ladderize}} function in the
#'   \code{\link[ape]{ape}} package performs a similar operation for objects
#'   of class \code{"phylo"}.
#' @examples
#'   newick <- "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
#'   dendro <- read.dendrogram(text = newick)
#'   plot(dendro, horiz = TRUE)
#'   dendro <- ladder(dendro, decreasing = TRUE)
#'   plot(dendro, horiz = TRUE)
#'
################################################################################
ladder <- function(x, decreasing = FALSE){
  ordernode <- function(tree, decreasing = FALSE){
    if(is.list(tree)){
      cladesizes <- sapply(tree, function(e) attr(e, "members"))
      if(any(cladesizes > 1)){
        cladeorder <- order(cladesizes, decreasing = decreasing)
        tree[] <- tree[cladeorder]
      }
    }
    return(tree)
  }
  ordernodes <- function(tree, decreasing = FALSE){
    if(is.list(tree)){
      tree[] <- lapply(tree, ordernode, decreasing = decreasing)
      tree[] <- lapply(tree, ordernodes, decreasing = decreasing)
    }
    return(tree)
  }
  x <- ordernodes(x, decreasing = decreasing)
  x <- ordernode(x, decreasing = decreasing)
  #shortcut in lieu of re-midpoint function *TODO*
  newick <- write.dendrogram(x)
  x <- read.dendrogram(text = newick)
  return(x)
}
