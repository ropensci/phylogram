#' Remove tree nodes by regular expression pattern matching.
#'
#' \code{"prune"} takes an object of class \code{"dendrogram"} and
#'   removes all branches whose branch labels match a given regular
#'   expression.
#'
#' @param tree an object of class \code{"dendrogram"}.
#' @param pattern a regular expression provided as a single string variable.
#' @param keep logical indicating whether the branches whose labels match
#'   the regexp pattern provided should be kept (TRUE) or discarded (FALSE;
#'   default)
#' @param untag logical (used only when keep = TRUE). Indicates whether
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
prune <- function(tree, pattern, keep = FALSE, untag = FALSE){
  collapse <- function(node, pattern, keep = FALSE){
    if(is.list(node)){
      childnames <- sapply(node, function(e) if(is.null(attr(e, "label"))) NA else attr(e, "label"))
      if(all(is.na(childnames))) return(node)
      condemned.leaves <- grepl(pattern, childnames)
      if(keep){
        condemned.leaves <- !condemned.leaves
        condemned.leaves[is.na(childnames)] <- FALSE
      }
      condemned.leaf <- match(TRUE, condemned.leaves)# just do one at a time so don't get singleton nodes
      if(!is.na(condemned.leaf)){
        tmpattr <- attributes(node) # cache internal node attributes
        node <- node[-condemned.leaf]
        tmpattr$members <- tmpattr$members - 1
        attributes(node) <- tmpattr
        if(length(node) == 1) node <- node[[1]] # delete inner node
      }
    }
    return(node)
  }
  prune1 <- function(tree, pattern, keep = FALSE){
    if(is.list(tree)){
      tree[] <- lapply(tree, collapse, pattern = pattern, keep = keep)
      tree[] <- lapply(tree, prune1, pattern = pattern, keep = keep)
    }
    return(tree)
  }
  fixmembers <- function(tree){
    attr(tree, "members") <- length(unlist(tree))
    return(tree)
  }
  #while(any(grepl(pattern, unlist(dendrapply(tree, function(e) attr(e, "label")))))){
  repeat{
    condemned.leaves <- grepl(pattern, unlist(dendrapply(tree, function(e) attr(e, "label"))))
    if(keep) condemned.leaves <- !condemned.leaves
    if(all(condemned.leaves)) return(structure(list(), class = "dendrogram", members = 0, height = 0))
    if(!any(condemned.leaves)) break
    tree <- prune1(tree, pattern = pattern, keep = keep)
    tree <- dendrapply(tree, fixmembers)
    tree <- collapse(tree, pattern = pattern, keep = keep)
  }
  newick <- write.dendrogram(tree)
  if(keep & untag) newick <- gsub(pattern, "", newick)
  tree <- read.dendrogram(text = newick)
  return(tree)
}



