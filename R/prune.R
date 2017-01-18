#' Remove tree branches based on regular expression matching.
#'
#' \code{"prune"} takes an object of class \code{"dendrogram"} and removes all branches whose
#' branch labels match a given regular expression.
#'
#' @param tree an object of class \code{"dendrogram"}.
#' @param pattern a regular expression provided as a single string variable.
#' @param keep logical indicating whether the branches whose labels match the regexp pattern
#'     provided should be kept (TRUE) or discarded (FALSE; default)
#' @param untag logical used only when keep = TRUE, indicating whether the pattern should be removed
#' from the branch labels.
#' @export
#'
prune <- function(tree, pattern, keep = FALSE, untag = FALSE){
  collapse <- function(tree, pattern, keep = FALSE){
    if(is.list(tree)){
      childnames <- sapply(tree, function(e) if(is.null(attr(e, "label"))) NA else attr(e, "label"))
      if(all(is.na(childnames))) return(tree)
      condemned.leaves <- grepl(pattern, childnames)
      if(keep){
        condemned.leaves <- !condemned.leaves
        condemned.leaves[is.na(childnames)] <- FALSE
      }
      condemned.leaf <- match(TRUE, condemned.leaves)# just do one at a time so don't get singleton trees
      if(!is.na(condemned.leaf)){
        tmpattr <- attributes(tree) # cache internal tree attributes
        tree <- tree[-condemned.leaf]
        tmpattr$members <- tmpattr$members - 1
        attributes(tree) <- tmpattr
        if(length(tree) == 1) tree <- tree[[1]] # delete inner node
      }
    }
    return(tree)
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
  #shortcut in lieu of re-midpoint function *TODO*
  newick <- write.dendrogram(tree)
  if(keep & untag) newick <- gsub(pattern, "", newick)
  tree <- read.dendrogram(text = newick)
  return(tree)
}



