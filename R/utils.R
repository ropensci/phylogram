#' Set dendrogram attributes for a nested list.
#'
#' \code{remidpoint} is a helper function used for manually creating
#'   \code{"dendrogram"} objects from nested lists. The function
#'   recursively assigns the necessary 'midpoint' and 'members'
#'   attributes at each node.
#'
#' @param x a nested list, possibly of class \code{"dendrogram"}
#' @param widths logical indicating whether the x coordinates should
#'   be included as "width" attributes at each node.
#' @return returns a nested list, or an object of class \code{"dendrogram"}
#'   depending on the class of the input object.
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
################################################################################
remidpoint <- function(x, widths = FALSE){
  isdendro <- inherits(x, "dendrogram")
  setnodeattr <- function(node, widths = FALSE){
    if(is.list(node)){
      cladesizes <- sapply(node, function(subnode){
        length(unlist(subnode, use.names = FALSE))
      })
      nclades <- length(cladesizes)
      if(widths){
        cladecounter <- 0
        for(i in 1:nclades){
          attr(node[[i]], "width") <- attr(node, "width") + cladecounter
          cladecounter <- cladecounter + cladesizes[i]
        }
      }
      attr(node, "members") <- sum(cladesizes)
      attr(node, "midpoint") <- ((cladesizes[1] - 1)/2 +
                                   #(cladesizes[1] +
                                   (sum(cladesizes[1:(nclades - 1)]) +
                                      (cladesizes[nclades] - 1)/2))/2
      if(widths) attr(node, "width") <- attr(node, "width") + attr(node, "midpoint")
      attr(node, "leaf") <- NULL
    }else{
      attr(node, "members") <- 1
      attr(node, "leaf") <- TRUE
    }
    return(node)
  }
  settreeattr <- function(tree, widths = FALSE){
    tree <- setnodeattr(tree, widths = widths)
    if(is.list(tree)) tree[] <- lapply(tree, settreeattr, widths = widths)
    return(tree)
  }
  if(widths) attr(x, "width") <- 1
  x <- settreeattr(x, widths = widths)
  #x <- setnodeattr(x)
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
#' @param shift either the character string "reset" (shift the graph so that
#'   the height of the farthest leaf from the root is zero), or a numeric value
#'   giving the amount to shift the graph along the primary axis.
#' @return returns an object of class \code{"dendrogram"}.
#' @author Shaun Wilkinson
#' @examples
#'   newick <- "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
#'   dendro <- read.dendrogram(text = newick)
#'   dendro <- reposition(dendro)
#'   plot(dendro, horiz = TRUE)
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
#' Make dendrogram ultrametric.
#'
#' This is a simple function that sets the 'height' attributes of
#'   all leaf nodes to  zero.
#'
#' @param x an object of class \code{"dendrogram"}.
#' @return returns an object of class \code{"dendrogram"}.
#' @author Shaun Wilkinson
#' @examples
#'   newick <- "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
#'   dendro <- read.dendrogram(text = newick)
#'   plot(dendro, horiz = TRUE)
#'   dendro <- ultrametricize(dendro)
#'   plot(dendro, horiz = TRUE)
################################################################################
ultrametricize <- function(x){
  if(!(inherits(x, "dendrogram"))) stop("x must be a 'dendrogram' object")
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
      if(all(sapply(node, function(e) !is.null(attr(e, "flag"))))){
        return(TRUE)
      }else return(FALSE)
    }else return(FALSE)
  }
  removeflags <- function(node){
    if(is.list(node)){
      for(i in seq_along(node)){
        attr(node[[i]], "flag") <- NULL
      }
    }
    return(node)
  }
  ultrametricize1 <- function(node){
    if(is.list(node)){
      if(checkflags(node)){
        childheights <- sapply(node, function(e) attr(e, "height"))
        attr(node, "height") <- max(childheights) + 1
        node <- removeflags(node)
        attr(node, "flag") <- TRUE
      }
    }
    return(node)
  }
  while(is.null(attr(x, "flag"))) x <- dendrapply(x, ultrametricize1)
  removeallflags <- function(node){
    attr(node, "flag") <- NULL
    return(node)
  }
  x <- dendrapply(x, removeallflags)
  return(x)
}
################################################################################
