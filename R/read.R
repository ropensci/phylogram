#' Read a dendrogram from parenthetic text.
#'
#' This function wraps the \code{\link[ape]{read.tree}} parser from the
#'   \code{\link[ape]{ape}} package to read a phylogenetic tree from
#'   parenthetic text in the Newick/New Hampshire format, and
#'   converts it to object of class "dendrogram".
#'
#' @param file character string giving a valid path to the file from
#'   which to read the data.
#' @param text optional character string in lieu of a "file" argument.
#'   If a text argument is provided instead of a file path, the data
#'   are read via a text connection.
#' @param ... further arguments to be passed to
#'   \code{\link[ape]{read.tree}} (which may then be passed on to
#'   \code{scan}).
#' @return an object of class \code{"dendrogram"}.
#' @author Shaun Wilkinson
#' @references
#'   Felsenstein J (1986) The Newick tree format.
#'   \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}
#'
#'   Olsen G (1990) Interpretation of the "Newick's 8:45" tree format standard.
#'   \url{http://evolution.genetics.washington.edu/phylip/newick_doc.html}
#'
#'   Paradis E, Claude J, Strimmer K, (2004) APE: analyses of phylogenetics
#'   and evolution in R language. \emph{Bioinformatics} \strong{20}, 289-290.
#'
#'   Paradis E (2008) Definition of Formats for Coding Phylogenetic Trees in R.
#'   \url{http://ape-package.ird.fr/misc/FormatTreeR_24Oct2012.pdf}
#'
#'   Paradis E (2012) Analysis of Phylogenetics and Evolution with R
#'   (Second Edition). Springer, New York.
#' @seealso
#'   \code{\link{write.dendrogram}} writes an object of
#'   class \code{"dendrogram"} to a Newick text string.
#'   The \code{\link[ape]{read.tree}} function in the
#'   \code{\link[ape]{ape}} package parses objects
#'   of class \code{"phylo"} and \code{"multiPhylo"}.
#' @examples
#'   x <- read.dendrogram(text = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
#'   plot(x, horiz = TRUE)
################################################################################
read.dendrogram <- function(file = "", text = NULL, ...){
  phy <- ape::read.tree(file = file, text = text, ... = ...)
  if(inherits(phy, "multiPhylo")) stop("Only single trees are supported\n")
  dnd <- as.dendrogram.phylo(phy)
  return(dnd)
}
################################################################################
