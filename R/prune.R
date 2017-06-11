#' Remove tree nodes by regular expression pattern matching.
#'
#' \code{"prune"} takes an object of class \code{"dendrogram"} and
#'   removes all branches whose branch labels match a given regular
#'   expression.
#'
#' @param tree an object of class \code{"dendrogram"}.
#' @param pattern a regular expression.
#' @param invert logical indicating whether the branches whose labels match
#'   the regular expression provided in "pattern" should be
#'   discarded and the others kept (FALSE; default) or vice versa.
#'   Nodes without "label" attributes are ignored.
#' @param untag logical (used only when invert = TRUE). Indicates whether
#'   the specified pattern should be removed from the branch labels in
#'   the returned object.
#' @param ... further arguments to be passed to \code{grepl} and \code{gsub}.
#' @return Returns an object of class \code{"dendrogram"}.
#'
#' @details This function recursively tests the "label" attribute of each
#'   dendrogram node (including non-leaf inner nodes if applicable) for
#'   the specified pattern, removing those that register a positive hit.
#'   Note that positive matching inner nodes are removed along with all of
#'   their sub-nodes, regardless of whether the "label" attributes of the
#'   sub-nodes match the pattern.
#'
#' @author Shaun Wilkinson
#'
#' @seealso The \code{\link[ape]{drop.tip}} function in the
#'   \code{\link[ape]{ape}} package performs a similar operation for objects
#'   of class \code{"phylo"}. See \code{\link{regex}} for help with
#'   compiling regular expressions.
#'
#' @examples
#'   x <- read.dendrogram(text = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
#'   plot(x, horiz = TRUE)
#'   x <- prune(x, pattern = "^A$")
#'   plot(x, horiz = TRUE)
#'
################################################################################
prune <- function(tree, pattern, invert = FALSE, untag = FALSE, ...){
  collapse <- function(node, pattern, invert = FALSE){
    if(is.list(node)){
      childnames <- sapply(node, function(e){
        if(is.null(attr(e, "label"))) NA else attr(e, "label")
      })
      if(all(is.na(childnames))) return(node)
      condemnedleaves <- grepl(pattern, childnames, ... = ...)
      if(invert){
        condemnedleaves <- !condemnedleaves
        condemnedleaves[is.na(childnames)] <- FALSE
      }
      condemnedleaf <- match(TRUE, condemnedleaves)
      # just do one at a time to prevent singleton nodes
      if(!is.na(condemnedleaf)){
        tmpattr <- attributes(node) # cache internal node attributes
        node <- node[-condemnedleaf]
        tmpattr$members <- tmpattr$members - 1
        attributes(node) <- tmpattr
        if(length(node) == 1) node <- node[[1]] # deletes inner node
      }
    }
    return(node)
  }
  collapser <- function(tree, pattern, invert = FALSE){
    tree <- collapse(tree, pattern = pattern, invert = invert)
    if(is.list(tree)) tree[] <- lapply(tree, collapser, pattern = pattern, invert = invert)
    return(tree)
  }
  fixmembers <- function(tree){
    attr(tree, "members") <- length(unlist(tree, use.names = FALSE))
    return(tree)
  }
  grepd <- function(pattern, x, invert = FALSE){
    # grepl the "label" attrs of a dendrogram
    temporaryflag <- FALSE
    grepnode <- function(node, pattern, invert = FALSE){
      if(!is.null(attr(node, "label"))){
        hit <- grepl(pattern, attr(node, "label"), ... = ...)
        if(invert) hit <- !hit
        if(hit) temporaryflag <<- TRUE
      }
      return(node)
    }
    tmp <- dendrapply(x, grepnode, pattern, invert = invert)
    out <- temporaryflag
    rm(temporaryflag)
    return(out)
  }
  repeat{
    condemnedleaves <- grepd(pattern, tree, invert = invert)
    if(!condemnedleaves) break
    if(length(tree) == 1){
      if(!is.null(attr(tree, "label"))){
        hit <- grepl(pattern, attr(tree, "label"), ... = ...)
        if(invert) hit <- !hit
        if(hit) return(structure(list(), class = "dendrogram", members = 0, height = 0))
      }
    }
    tree <- collapser(tree, pattern = pattern, invert = invert)
    tree <- dendrapply(tree, fixmembers)
  }
  if(length(tree) == 0) return(NULL)
  if(length(tree) == 1 & is.list(tree)) tree <- tree[[1]]
  if(invert & untag){
    untagnode <- function(node, tag){
      if(!is.null(attr(node, "label"))){
        attr(node, "label") <- gsub(tag, "", attr(node, "label"), ... = ...)
      }
      return(node)
    }
    tree <- dendrapply(tree, untagnode, tag = pattern)
  }
  tree <- remidpoint(tree)
  return(tree)
}
################################################################################
