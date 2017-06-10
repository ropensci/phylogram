#' Reorder tree branches in ladderized pattern.
#'
#' This function ladderizes the branches of a \code{dendrogram} object
#'   to aid in visual interpretation.
#'
#' @param x an object of class \code{"dendrogram"}.
#' @param decreasing logical indicating whether the tree should be
#'   ladderized upwards or downwards. Defaults to FALSE (downwards).
#' @return Returns an object of class \code{dendrogram}.
#' @details This function is the \code{dendrogram} analogue of the
#'   \code{\link[ape]{ladderize}} function in the \code{\link[ape]{ape}}
#'   package (Paradis et al 2004, 2012).
#' @author Shaun Wilkinson
#' @references
#'   Paradis E, Claude J, Strimmer K, (2004) APE: analyses of phylogenetics
#'   and evolution in R language. \emph{Bioinformatics} \strong{20}, 289-290.
#'
#'   Paradis E (2012) Analysis of Phylogenetics and Evolution with R
#'   (Second Edition). Springer, New York.
#'
#' @seealso The \code{\link[ape]{ladderize}} function in the
#'   \code{\link[ape]{ape}} package performs a similar operation for objects
#'   of class \code{"phylo"}.
#'
#' @examples
#'   x <- read.dendrogram(text = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
#'   plot(x, horiz = TRUE)
#'   x <- ladder(x, decreasing = TRUE)
#'   plot(x, horiz = TRUE)
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
  x <- remidpoint(x)
  return(x)
}
################################################################################
