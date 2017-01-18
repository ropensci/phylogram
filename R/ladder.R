#' Reorder the branches of a dendrogram in a ladderized fashion.
#'
#' \code{"ladder"} takes an object of class "dendrogram" and ladderizes
#' the branches for easier visual interpretation.
#'
#' @param x an object of class \code{"dendrogram"}.
#' @param decreasing logical indicating whether the tree should be ladderized upwards
#' (FALSE; default) or downwards (TRUE).
#' @export
#'
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
