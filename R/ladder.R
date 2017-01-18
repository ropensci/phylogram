ladder <- function(tree, decreasing = FALSE){
  ordernode <- function(tree, decreasing = FALSE){
    if(is.list(tree)){
      cladesizes <- sapply(tree, function(e) attr(e, "members"))
      if(any(cladesizes > 1)){
        cladeorder <- order(cladesizes, decreasing = decreasing)
        tree[] <- tree[cladeorder]
      }
    }
    return(tree)
  }
  ordernodes <- function(tree, decreasing = FALSE){
    if(is.list(tree)){
      tree[] <- lapply(tree, ordernode, decreasing = decreasing)
      tree[] <- lapply(tree, ordernodes, decreasing = decreasing)
    }
    return(tree)
  }
  tree <- ordernodes(tree, decreasing = decreasing)
  tree <- ordernode(tree, decreasing = decreasing)
  #shortcut in lieu of re-midpoint function *TODO*
  newick <- write.dendrogram(tree)
  tree <- read.dendrogram(text = newick)
  return(tree)
}
