#' Remove tree nodes with regular expression pattern matching.
#'
#' \code{"prune"} takes an object of class \code{"dendrogram"} and
#'   removes all branches whose branch labels match a given regular
#'   expression.
#'
#' @param tree an object of class \code{"dendrogram"}.
#' @param pattern a regular expression provided as a single string variable.
#' @param invert logical indicating whether the branches whose labels match
#'   the regular expression provided in "pattern" should be
#'   discarded (FALSE; default) or kept.
#'   Nodes with no "label" attributes are ignored.
#' @param untag logical (used only when invert = TRUE). Indicates whether
#'   the specified pattern should be removed from the branch labels in
#'   the returned object.
#' @return an object of class \code{"dendrogram"}.
#' @details TBA.
#' @author Shaun Wilkinson
#' @seealso The \code{\link[ape]{drop.tip}} function in the
#'   \code{\link[ape]{ape}} package performs a similar operation for objects
#'   of class \code{"phylo"}.
#' @examples
#'   newick <- "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
#'   dendro <- read.dendrogram(text = newick)
#'   plot(dendro, horiz = TRUE)
#'   dendro <- prune(dendro, pattern = "^A$")
#'   plot(dendro, horiz = TRUE)
#'
################################################################################
prune <- function(tree, pattern, invert = FALSE, untag = FALSE){
  collapse <- function(node, pattern, invert = FALSE){
    if(is.list(node)){
      childnames <- sapply(node, function(e){
        if(is.null(attr(e, "label"))) NA else attr(e, "label")
      })
      if(all(is.na(childnames))) return(node)
      condemnedleaves <- grepl(pattern, childnames)
      if(invert){
        condemnedleaves <- !condemnedleaves
        condemnedleaves[is.na(childnames)] <- FALSE
      }
      condemnedleaf <- match(TRUE, condemnedleaves)
      # just do one at a time so not to get singleton nodes
      if(!is.na(condemnedleaf)){
        tmpattr <- attributes(node) # cache internal node attributes
        node <- node[-condemnedleaf]
        tmpattr$members <- tmpattr$members - 1 # necessary?
        attributes(node) <- tmpattr
        if(length(node) == 1) node <- node[[1]] # deletes inner node
      }
    }
    return(node)
  }
  prune1 <- function(tree, pattern, invert = FALSE){
    if(is.list(tree)){
      tree[] <- lapply(tree, collapse, pattern = pattern, invert = invert)
      tree[] <- lapply(tree, prune1, pattern = pattern, invert = invert)
    }
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
        hit <- grepl(pattern, attr(node, "label"))
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
    # condemnedleaves <- grepl(pattern, unlist(dendrapply(tree, function(e) attr(e, "label"))))
    # if(invert) condemnedleaves <- !condemnedleaves
    # if(all(condemnedleaves)) {
    #   return(structure(list(), class = "dendrogram", members = 0, height = 0))
    # }
    # if(!any(condemnedleaves)) break
    if(!condemnedleaves) break
    if(length(tree) == 1){
      if(!is.null(attr(tree, "label"))){
        hit <- grepl(pattern, attr(tree, "label"))
        if(invert) hit <- !hit
        if(hit) return(structure(list(), class = "dendrogram", members = 0, height = 0))
      }
    }
    tree <- prune1(tree, pattern = pattern, invert = invert)
    tree <- dendrapply(tree, fixmembers)
    tree <- collapse(tree, pattern = pattern, invert = invert)
  }
  if(length(tree) == 0) return(NULL)
  newick <- write.dendrogram(tree)
  if(invert & untag) newick <- gsub(pattern, "", newick)
  tree <- read.dendrogram(text = newick)
  return(tree)
}



