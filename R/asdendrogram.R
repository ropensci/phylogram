#' Convert a "phylo" object to a dendrogram.
#'
#' This function converts a "phylo" object (Paradis et al 2004)
#'   to a dendrogram.
#'
#' @param object an object of class "phylo".
#' @param ... further arguments to be passed between methods.
#' @return an object of class "dendrogram".
#' @author Shaun Wilkinson
#' @references
#'   Paradis E, Claude J, Strimmer K, (2004) APE: analyses of phylogenetics
#'   and evolution in R language. \emph{Bioinformatics} \strong{20}, 289-290.
#'
#'   Paradis E (2008) Definition of Formats for Coding Phylogenetic Trees in R.
#'   \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}
#'
#'   Paradis E (2012) Analysis of Phylogenetics and Evolution with R
#'   (Second Edition). Springer, New York.
#' @examples
#'   newick <- "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
#'   x <- read.dendrogram(text = newick)
#'   y <- as.phylo(x)
#'   z <- as.dendrogram(y)
#'   identical(x, z)
#' @method as.dendrogram phylo
################################################################################
as.dendrogram.phylo <- function(object, ...){
  dots <- list(...)
  hasedges <- !is.null(object$edge.length)
  if(!hasedges) object$edge.length <- rep(1, nrow(object$edge))
  hasnodelabs <- !is.null(object$node.label)
  if(hasnodelabs) stopifnot(length(object$node.label) == object$Nnode)
  minheight <- 1E07
  tree <- 1L
  attr(tree, "height") <- 0
  attr(tree, "node") <- length(object$tip.label) + 1L
  lookups <- lapply(seq(1, max(object$edge[, 1])), function(e) integer(0))
  for(i in seq_len(nrow(object$edge))){
    lookups[[object$edge[i, 1]]] <- c(lookups[[object$edge[i, 1]]], i)
  }
  splitnode <- function(node, object){
    nodenumber <- attr(node, "node")
    if(is.null(nodenumber)) return(node) # for tips only
    nodeheight <- attr(node, "height")
    if(nodenumber > length(object$tip.label)){ #inner nodes
      indices <- lookups[[nodenumber]]
      childnodes <- object$edge[indices, 2]
      nchildren <- length(childnodes)
      childheights <- nodeheight - object$edge.length[indices]
      node <- vector(mode = "list", length = nchildren)
      attr(node, "node") <- nodenumber
      attr(node, "height") <- nodeheight
      for(i in seq_along(node)){
        node[[i]] <- 1L
        attr(node[[i]], "height") <- childheights[i]
        attr(node[[i]], "node") <- childnodes[i]
      }
    }else{ #leaf nodes
      node <- nodenumber
      if(nodeheight < minheight) minheight <<- nodeheight
      attr(node, "height") <- nodeheight
      attr(node, "label") <- object$tip.label[nodenumber]
      attr(node, "leaf") <- TRUE
    }
    return(node)
  }
  splittree <- function(tree, object){
    tree <- splitnode(tree, object)
    if(is.list(tree)) tree[] <- lapply(tree, splittree, object)
    return(tree)
  }
  tree <- splittree(tree, object)
  tree <- remidpoint(tree)
  class(tree) <- "dendrogram"
  tree <- if(hasedges){
    reposition(tree, shift = -1 * minheight)
  }else{
    as.cladogram(tree)
  }
  if(hasnodelabs){
    inc <- 1
    fun <- function(node, object){
      if(is.list(node)){
        attr(node, "label") <- object$node.label[inc]
        inc <<- inc + 1
      }
      return(node)
    }
    tree <- dendrapply(tree, fun, object)
  }
  if(!is.null(object$root.edge)) attr(tree, "edge.root") <- object$root.edge
  return(tree)
}
################################################################################
