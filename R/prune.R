#' Remove tree nodes by regular expression pattern matching.
#'
#' \code{"prune"} takes an object of class \code{"dendrogram"} and
#'   removes all branches whose branch labels match a given regular
#'   expression.
#'
#' @param tree an object of class \code{"dendrogram"}.
#' @param pattern a regular expression.
#' @param keep logical indicating whether the nodes whose labels match
#'   the regular expression provided in "pattern" should be
#'   kept and the remainder discarded. Defaults to FALSE.
#'   Note that nodes without "label" attributes are ignored.
#' @param ... further arguments to be passed to \code{grepl} and \code{gsub}.
#' @return Returns an object of class \code{"dendrogram"}.
#' @details This function recursively tests the "label" attribute of each
#'   dendrogram node (including non-leaf inner nodes if applicable) for
#'   the specified pattern, removing those that register a positive hit.
#'   Note that positive matching inner nodes are removed along with all of
#'   their sub-nodes, regardless of whether the "label" attributes of the
#'   sub-nodes match the pattern.
#' @author Shaun Wilkinson
#' @seealso The \code{\link[ape]{drop.tip}} function in the
#'   \code{\link[ape]{ape}} package performs a similar operation for objects
#'   of class \code{"phylo"}. See \code{\link{regex}} for help with
#'   compiling regular expressions.
#' @examples
#'   x <- read.dendrogram(text = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
#'   plot(x, horiz = TRUE)
#'   x <- prune(x, pattern = "^A$")
#'   plot(x, horiz = TRUE)
################################################################################
prune <- function(tree, pattern, keep = FALSE, ...){
  if(!is.null(attr(tree, "label"))){
    if(!keep & grepl(pattern, attr(tree, "label"), ... = ...)){
      return(structure(list(), class = "dendrogram", members = 0, height = 0))
    }
  }
  # if(length(pattern) > 1) pattern <- paste0(pattern, collapse = "|")
  collapse <- function(node, pattern, keep = FALSE){
    if(is.list(node)){
      lab <- function(e) if(is.null(attr(e, "label"))) "" else attr(e, "label")
      childnames <- vapply(node, lab, "")
      if(all(childnames == "")) return(node)
      condemnedleaves <- grepl(pattern, childnames, ... = ...)
      if(keep){
        condemnedleaves <- !condemnedleaves
        condemnedleaves[childnames == ""] <- FALSE
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
  collapser <- function(tree, pattern, keep = FALSE){
    tree <- collapse(tree, pattern = pattern, keep = keep)
    if(is.list(tree)){
      tree[] <- lapply(tree, collapser, pattern = pattern, keep = keep)
    }
    return(tree)
  }
  fixmembers <- function(tree){
    attr(tree, "members") <- length(unlist(tree, use.names = FALSE))
    return(tree)
  }
  grepd <- function(pattern, x, keep = FALSE){
    # grepl the "label" attrs of a dendrogram
    out <- FALSE
    grepnode <- function(node, pattern, keep = FALSE, ...){
      if(!is.null(attr(node, "label"))){
        hit <- grepl(pattern, attr(node, "label"), ... = ...)
        if(keep) hit <- !hit
        if(hit) out <<- TRUE
      }
      return(node)
    }
    tmp <- dendrapply(x, grepnode, pattern, keep = keep, ... = ...)
    return(out)
  }
  repeat{
    condemnedleaves <- grepd(pattern, tree, keep = keep)
    if(!condemnedleaves) break
    if(length(tree) == 1){
      if(!is.null(attr(tree, "label"))){
        hit <- grepl(pattern, attr(tree, "label"), ... = ...)
        if(keep) hit <- !hit
        if(hit){
          return(structure(list(), class = "dendrogram",
                           members = 0, height = 0))
        }
      }
    }
    tree <- collapser(tree, pattern = pattern, keep = keep)
    tree <- dendrapply(tree, fixmembers)
  }
  if(length(tree) == 0) return(NULL)
  if(length(tree) == 1 & is.list(tree)) tree <- tree[[1]]
  tree <- remidpoint(tree)
  return(tree)
}
################################################################################
