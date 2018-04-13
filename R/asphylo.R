#' Convert dendrogram to "phylo" object.
#'
#' Functions for converting dendrograms to "phylo" objects
#'   and \emph{vice versa}.
#'
#' @param x a dendrogram object.
#' @param ... further arguments to be passed between methods.
#' @return an object of class "phylo".
#' @details These functions currently work by temporarily writing a tree to Newick
#'   text and then parsing the string using either \code{\link{read.dendrogram}}
#'   or \code{\link[ape]{read.tree}}. A faster implementation that avoids the
#'   transformation to text strings will be available in a future version.
#' @author Shaun Wilkinson
#' @examples
#'   newick <- "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
#'   x <- read.dendrogram(text = newick)
#'   y <- as.phylo(x)
#'   z <- as.dendrogram(y)
#'   identical(x, z)
#' @method as.phylo dendrogram
################################################################################
as.phylo.dendrogram <- function(x, ...){
  ape::read.tree(text = write.dendrogram(x), ... = ...)
}
################################################################################
