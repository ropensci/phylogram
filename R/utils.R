#' Set dendrogram attributes for a nested list.
#'
#' \code{remidpoint} is a helper function used for manually creating
#'   \code{"dendrogram"} objects from nested lists. The function
#'   recursively assigns the necessary 'midpoint', 'members', and
#'   'leaf' attributes at each node.
#'
#' @param x a nested list, possibly of class \code{"dendrogram"}
#' @return returns a nested list, or an object of class \code{"dendrogram"}
#'   depending on the class of the input object.
#' @author Shaun Wilkinson
#' @examples
#'   ## manually create a small dendrogram with three members, A, B and C
#'   x <- list("A", list("B", "C"))
#'   attr(x[[1]], "leaf") <- TRUE
#'   attr(x[[2]][[1]], "leaf") <- TRUE
#'   attr(x[[2]][[2]], "leaf") <- TRUE
#'   attr(x[[1]], "label") <- "A"
#'   attr(x[[2]][[1]], "label") <- "B"
#'   attr(x[[2]][[2]], "label") <- "C"
#'   attr(x, "height") <- 1
#'   attr(x[[1]], "height") <- 0
#'   attr(x[[2]], "height") <- 0.5
#'   attr(x[[2]][[1]], "height") <- 0
#'   attr(x[[2]][[2]], "height") <- 0
#'   x <- remidpoint(x)
#'   class(x) <- "dendrogram"
#'   plot(x, horiz = TRUE)
################################################################################
remidpoint <- function(x){
  isdendro <- inherits(x, "dendrogram")
  setnodeattr <- function(node){
    if(is.list(node)){
      get_cladesizes <- function(z) length(unlist(z, use.names = FALSE))
      cladesizes <- vapply(node, get_cladesizes, 0)
      nclades <- length(cladesizes)
      attr(node, "members") <- sum(cladesizes)
      attr(node, "midpoint") <- ((cladesizes[1] - 1)/2 +
                                   #(cladesizes[1] +
                                   (sum(cladesizes[1:(nclades - 1)]) +
                                      (cladesizes[nclades] - 1)/2))/2
      attr(node, "leaf") <- NULL
    }else{
      attr(node, "members") <- 1
      attr(node, "leaf") <- TRUE
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
#' Reset dendrogram height attributes.
#'
#' \code{reposition} is a helper function used for manually creating
#'   \code{"dendrogram"} objects from nested lists. The function
#'   recursively reassigns the 'height' attributes at each node by
#'   a given constant.
#'
#' @param x an object of class \code{"dendrogram"}.
#' @param shift either the character string "reset" (shift the graph so that
#'   the height of the farthest leaf from the root is zero), or a numeric value
#'   giving the amount to shift the graph along the primary axis.
#' @return Returns an object of class \code{"dendrogram"}.
#' @author Shaun Wilkinson
#' @examples
#'   x <- read.dendrogram(text = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
#'   plot(x, horiz = TRUE)
#'   x <- reposition(x)
#'   plot(x, horiz = TRUE)
################################################################################
reposition <- function(x, shift = "reset"){
  if(!(inherits(x, "dendrogram"))) stop("Input object must be of class
                                        'dendrogram'")
  if(identical(shift, "reset")){
    shift <- -1 * min(unlist(dendrapply(x, attr, "height"), use.names = FALSE))
  }
  reposition1 <- function(node, shift){ # node is a dendrogram
    attr(node, "height") <- attr(node, "height") + shift
    return(node)
  }
  x <- dendrapply(x, reposition1, shift = shift)
  return(x)
}
################################################################################
#' Apply unweighted branch lengths.
#'
#' This function sets the 'height' attributes of all leaf nodes to zero and
#'   progressively resets the heights of the inner nodes by single incremental
#'   units in a bottom-up fashion.
#'
#' @param x an object of class \code{"dendrogram"}.
#' @return an object of class \code{"dendrogram"}.
#' @author Shaun Wilkinson
#' @examples
#'   x <- read.dendrogram(text = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
#'   plot(x, horiz = TRUE)
#'   x <- as.cladogram(x)
#'   plot(x, horiz = TRUE)
################################################################################
as.cladogram <- function(x){
  if(!(inherits(x, "dendrogram"))){
    x <- remidpoint(x)
    class(x) <- "dendrogram"
  } # stop("x must be a 'dendrogram' object")
  placeflags <- function(node){
    if(is.leaf(node)){
      attr(node, "height") <- 0
      attr(node, "flag") <- TRUE
    }
    return(node)
  }
  x <- dendrapply(x, placeflags)
  checkflags <- function(node){
    if(is.list(node)){
      has_flag <- function(e) !is.null(attr(e, "flag"))
      return(all(vapply(node, has_flag, logical(1))))
    }else{
      return(FALSE)
    }
  }
  removeflags <- function(node){
    if(is.list(node)){
      for(i in seq_along(node)){
        attr(node[[i]], "flag") <- NULL
      }
    }
    return(node)
  }
  ac1 <- function(node){
    if(is.list(node)){
      if(checkflags(node)){
        childheights <- vapply(node, attr, 0, "height")
        attr(node, "height") <- max(childheights) + 1
        node <- removeflags(node)
        attr(node, "flag") <- TRUE
      }
    }
    return(node)
  }
  while(is.null(attr(x, "flag"))) x <- dendrapply(x, ac1)
  removeallflags <- function(node){
    attr(node, "flag") <- NULL
    return(node)
  }
  x <- dendrapply(x, removeallflags)
  return(x)
}
################################################################################
