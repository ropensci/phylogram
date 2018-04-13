#' Convert "phylo" objects to dendrograms.
#'
#' These functions are used for converting dendrograms to "phylo" objects
#'   and \emph{vice versa}.
#'
#' @param object a "phylo" object.
#' @param ... further arguments to be passed between methods.
#' @return an object of class "dendrogram".
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
#' @method as.dendrogram phylo
################################################################################
as.dendrogram.phylo <- function(object, ...){
  read.dendrogram(text = ape::write.tree(object), ... = ...)
}
################################################################################
