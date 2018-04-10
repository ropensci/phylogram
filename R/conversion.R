#' Conversion between "dendrogram" and "phylo" objects.
#'
#' These functions are used for converting dendrograms to "phylo" objects
#'   and \emph{vice versa}.
#'
#' @param x a dendrogram or "phylo" object
#' @return an object of class "dendrogram" (for \code{as.dendrogram.phylo}) or
#'   "phylo" (for \code{as.phylo.dendrogram}).
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
#' @name conversion
################################################################################
as.phylo <- function(x) UseMethod("as.phylo")
################################################################################
#' @rdname conversion
################################################################################
as.dendrogram <- function(x) UseMethod("as.dendrogram")
################################################################################
#' @rdname conversion
################################################################################
as.dendrogram.phylo <- function(x) read.dendrogram(text=ape::write.tree(x))
################################################################################
#' @rdname conversion
################################################################################
as.phylo.dendrogram <- function(x) ape::read.tree(text=write.dendrogram(x))
################################################################################
