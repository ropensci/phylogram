#' Convert a dendrogram to a "phylo" object.
#'
#' This function converts a dendrogram into an object of class "phylo"
#'   (see Paradis et al 2004).
#'
#' @param x a dendrogram.
#' @param ... further arguments to be passed between methods.
#' @return an object of class "phylo".
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
#' @method as.phylo dendrogram
################################################################################
as.phylo.dendrogram <- function(x, ...){
  dots <- list(...)
  nmb <- attr(x, "members")
  edgemat <- matrix(0L, nrow = 2 * nmb, ncol = 2)
  edgelens <- numeric(2 * nmb)
  tiplabs <- character(nmb)
  inc <- nmb + 1L # inner node counter
  lnc <- 1L # leaf node counter
  ec <- 0L # edge counter
  ## attach special attributes to each node
  fun <- function(node){
    if(is.list(node)){
      for(i in seq_along(node)){
        attr(node, "edgelengths") <- attr(node, "height") -
          vapply(node, attr, 0, "height")
      }
      attr(node, "node") <- inc
      inc <<- inc + 1L
    }else{
      attr(node, "node") <- lnc
      lnc <<- lnc + 1L
    }
    attr(node, "order") <- ec
    ec <<- ec + 1L
    return(node)
  }
  x <- dendrapply(x, fun)
  ri <- integer(ec - 1) # reordering indices
  ec <- 1 # edge counter
  lnc <- 1 # leaf node counter
  ## build edge matrix etc
  fun <- function(node){
    if(is.list(node)){
      for(i in seq_along(node)){
        edgemat[ec, 1] <<- attr(node, "node")
        edgemat[ec, 2] <<- attr(node[[i]], "node")
        edgelens[ec] <<- attr(node, "edgelengths")[i]
        ri[ec] <<- attr(node[[i]], "order")
        ec <<- ec + 1
      }
    }else{
      tiplabs[lnc] <<- attr(node, "label")
      lnc <<- lnc + 1
    }
    return(node)
  }
  x <- dendrapply(x, fun)
  ri <- order(ri)
  keeps <- c(!logical(ec - 1L), logical(length(edgelens) - ec + 1L))
  edgemat <- edgemat[keeps, ]
  edgemat <- edgemat[ri, ]
  edgelens <- edgelens[keeps]
  edgelens <- edgelens[ri]
  nnodes <- inc - lnc
  phy <- list(edge = edgemat, edge.length = edgelens, Nnode = nnodes,
              tip.label = tiplabs)
  class(phy) <- "phylo"
  attr(phy, "order") <- "cladewise"
  hasnodelabs <- is.list(x) & !is.null(attr(x, "label"))
  if(hasnodelabs){
    inc <- 1
    nodelabs <- character(nmb)
    fun <- function(node){
      if(is.list(node)){
        if(is.null(attr(node, "label"))) attr(node, "label") <- ""
        nodelabs[inc] <<- attr(node, "label")
        inc <<- inc + 1
      }
      return(node)
    }
    x <- dendrapply(x, fun)
    phy$node.label <- nodelabs[1:nnodes]
  }
  if(!is.null(attr(x, "edge.root"))){
    phy$root.edge <- attr(x, "edge.root")
  }
  return(phy)
}
################################################################################
