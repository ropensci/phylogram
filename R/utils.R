#' Set dendrogram attributes for a nested list.
#'
#' \code{remidpoint} is a helper function used for manually creating
#'   \code{"dendrogram"} objects from nested lists. The function
#'   recursively assigns the necessary 'midpoint' and 'members'
#'   attributes at each node.
#'
#' @param x a nested list, possibly of class \code{"dendrogram"}
#' @return returns a nested list, or an object of class \code{"dendrogram"}
#'   depending on the class of the input object.
#' @details TBA.
#' @author Shaun Wilkinson
#' @examples
#'   futuredendro <- list("A", list("B", "C"))
#'   attr(futuredendro[[1]], "leaf") <- TRUE
#'   attr(futuredendro[[2]][[1]], "leaf") <- TRUE
#'   attr(futuredendro[[2]][[2]], "leaf") <- TRUE
#'   attr(futuredendro[[1]], "label") <- "A"
#'   attr(futuredendro[[2]][[1]], "label") <- "B"
#'   attr(futuredendro[[2]][[2]], "label") <- "C"
#'   attr(futuredendro, "height") <- 1
#'   attr(futuredendro[[1]], "height") <- 0
#'   attr(futuredendro[[2]], "height") <- 0.5
#'   attr(futuredendro[[2]][[1]], "height") <- 0
#'   attr(futuredendro[[2]][[2]], "height") <- 0
#'   dendro <- remidpoint(futuredendro)
#'   class(dendro) <- "dendrogram"
#'   plot(dendro, horiz = TRUE)
#'
################################################################################
remidpoint <- function(x){
  isdendro <- inherits(x, "dendrogram")
  setnodeattr <- function(node){
    if(is.list(node)){
      cladesizes <- sapply(node, function(subnode) length(unlist(subnode)))
      attr(node, "members") <- sum(cladesizes)
      attr(node, "midpoint") <- ((cladesizes[1] - 1)/2 +
                                   (cladesizes[1] +
                                      (cladesizes[2] - 1)/2))/2
    }else{
      attr(node, "members") <- 1
    }
    return(node)
  }
  settreeattr <- function(tree){
    tree <- setnodeattr(tree)
    if(is.list(tree)) tree[] <- lapply(tree, settreeattr)
    return(tree)
  }
  x <- settreeattr(x)
  if(isdendro) class(x) <- "dendrogram"
  return(x)
}
################################################################################

#' Reset dendrogram height to zero.
#'
#' \code{reposition} is a helper function used for manually creating
#'   \code{"dendrogram"} objects from nested lists. The function
#'   recursively assigns the necessary 'height' attributes at each node
#'   so that the height of the furthest leaf from the root node is zero.
#'
#' @param x an object of class \code{"dendrogram"}.
#' @return returns an object of class \code{"dendrogram"}.
#' @details TBA.
#' @author Shaun Wilkinson
#' @examples
#'   newick <- "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
#'   dendro <- read.dendrogram(text = newick)
#'   dendro <- reposition(dendro)
#'   plot(dendro, horiz = TRUE)
#'
################################################################################
reposition <- function(x){
  if(!(inherits(x, "dendrogram"))) stop("Input object must be of class
                                        'dendrogram'")
  minhgt <- min(unlist(dendrapply(x, attr, "height")))
  reposition1 <- function(node, minhgt){ # node is a dendrogram
    attr(node, "height") <- attr(node, "height") - minhgt
    return(node)
  }
  x <- dendrapply(x, reposition1, minhgt = minhgt)
  return(x)
}
################################################################################

#' Make dendrogram ultrametric.
#'
#' This is a simple function that sets the 'height' attributes of
#'   all leaf nodes to  zero.
#'
#' @param x an object of class \code{"dendrogram"}.
#' @return returns an object of class \code{"dendrogram"}.
#' @details TBA.
#' @author Shaun Wilkinson
#' @examples
#'   newick <- "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
#'   dendro <- read.dendrogram(text = newick)
#'   plot(dendro, horiz = TRUE)
#'   dendro <- ultrametricize(dendro)
#'   plot(dendro, horiz = TRUE)
#'
################################################################################
ultrametricize <- function(x){
  if(!(inherits(x, "dendrogram"))) stop("Input object must be of class
                                        'dendrogram'")
  ultrametricize1 <- function(node){
    if(is.leaf(node)) attr(node, "height") <- 0
    return(node)
  }
  x <- dendrapply(x, ultrametricize1)
  return(x)
}
################################################################################
