#' Internal phylogram functions
#'
remidpoint <- function(x){
  # a dendrogram whose nodes need re-midpointing
  setnodeattr <- function(node){
    if(is.list(node)){
      cladesizes <- sapply(node, function(subnode) length(unlist(subnode)))
      attr(node, "members") <- sum(cladesizes)
      attr(node, "midpoint") <- ((cladesizes[1] - 1)/2 + (cladesizes[1] + (cladesizes[2] - 1)/2))/2
    }else{
      attr(node, "members") <- 1
    }
    return(node)
  }
  settreeattr <- function(tree){
    tree <- setnodeattr(tree)
    if(is.list(tree)) tree[] <- lapply(tree, settreeattr)
    return(tree)
  }
  x <- settreeattr(x)
  return(x)
}

reposition <- function(x){
  minhgt <- min(unlist(dendrapply(x, attr, "height")))
  reposition1 <- function(node, minhgt){ # node is a dendrogram
    attr(node, "height") <- attr(node, "height") - minhgt
    return(node)
  }
  x <- dendrapply(x, reposition1, minhgt = minhgt)
  return(x)
}

ultrametricize <- function(x){
  ultrametricize1 <- function(node){
    if(is.leaf(node)) attr(node, "height") <- 0
    return(node)
  }
  x <- dendrapply(x, ultrametricize)
  return(x)
}

