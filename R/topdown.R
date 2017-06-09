#' Top down phylogenetic trees.
#'
#' Create phylogenetic trees by successively splitting the sequence dataset
#'   into smaller and smaller subsets.
#'
#' @param x a list or matrix of sequences, possibly an object of class
#'   \code{"DNAbin"} or \code{"AAbin"}. Alternatively can be an object of class
#'   \code{kcounts}.
#' @param k integer representing the the k-mer size to be used for calculating
#'   the distance matrix.
#' @param residues either NULL (default; emitted residues are automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Specifying the residue alphabet is therefore
#'   recommended unless the sequence list is a "DNAbin" or "AAbin" object.
#' @param gap the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} or \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @param ... further arguments to be passed to \code{kmeans} (not including
#'   \code{centers}).
#' @return an object of class \code{"dendrogram"}.
#' @details This function creates phylogenetic trees by successively splitting
#'   the sequence dataset into smaller and smaller subsets (a.k.a. recursive
#'   partitioning). This is a "top down" approach to tree-building, as
#'   opposed to agglomerative "bottom up" methods such as neighbour joining
#'   and UPGMA. It is useful for large sequence sets (> 20, 000) since the
#'   usual N x N distance matrix is not required (where N is the number of
#'   sequences). Instead, an N x log(N, 2)^2
#'   distance matrix is first derived from the input sequences using the
#'   \code{\link{mbed}} function of Blacksheilds et al. 2010 and the k-mer
#'   distance measure of (Edgar 2004). This is done by first selecting
#'   log(N, 2)^2 'seed' sequences (selected randomly by default)
#'   and calculating the distance of each sequence to each of the seeds.
#'   This requires substantially less memory and time than N x N matrix
#'   computation. The 'embedded' sequences
#'   are then split recursively using the k-means (clusters = 2) technique.
#' @author Shaun Wilkinson
#' @references
#'   Blackshields G, Sievers F, Shi W, Wilm A, Higgins DG (2010) Sequence embedding
#'   for fast construction of guide trees for multiple sequence alignment.
#'   \emph{Algorithms for Molecular Biology}, \strong{5}, 21.
#'
#'   Edgar RC (2004) Local homology recognition and distance measures in
#'   linear time using compressed amino acid alphabets.
#'   \emph{Nucleic Acids Research}, \strong{32}, 380-385.
#' @examples
#' \dontrun{
#' ## Cluster the woodmouse dataset from the ape package
#' library(ape)
#' data(woodmouse)
#' ## trim gappy ends for global alignment
#' woodmouse <- woodmouse[, apply(woodmouse, 2, function(v) !any(v == 0xf0))]
#' ## build tree divisively
#' set.seed(999)
#' woodmouse.tree <- topdown(woodmouse)
#' ## plot tree
#' op <- par(no.readonly = TRUE)
#' par(mar = c(5, 2, 4, 8) + 0.1)
#' plot(woodmouse.tree, main = "Woodmouse phylogeny", horiz = TRUE)
#' par(op)
#' }
################################################################################
topdown <- function(x, k = 5, residues = NULL, gap = "-", ...){
  # seeds = NULL, weighted = TRUE,
  # first embed the seqs in a N x log(N, 2)^2 distmat as in Blackshields 2010
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  residues <- .alphadetect(x, residues = residues, gap = gap)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  if(is.matrix(x)) x <- .unalign(x, gap = gap)
  nseq <- length(x)
  if(nseq == 1) stop("Only a single sequence provided")
    # M <- mbed(x, seeds = seeds, k = k, residues = residues, gap = gap,
    #           counts = weighted) # kcounts only required for weighted option
  if(is.null(names(x))) names(x) <- paste0("S", 1:nseq)
  catchnames <- names(x)
  hashes <- .digest(x, simplify = TRUE)
  duplicates <- duplicated(hashes)
  nuseq <- sum(!duplicates)
  if(any(duplicates)){
    pointers <- integer(length(x))
    dupehashes <- hashes[duplicates]
    uniquehashes <- hashes[!duplicates]
    pointers[!duplicates] <- seq_along(uniquehashes)
    pd <- integer(length(dupehashes))
    for(i in unique(dupehashes)) pd[dupehashes == i] <- match(i, uniquehashes)
    pointers[duplicates] <- pd
    x <- x[!duplicates]
  }else{
    pointers <- seq_along(x)
  }
  kcounts <- kcount(x, k = k, residues = residues, gap = gap)
  #duplicates <- attr(M, "duplicates")
  # pointers <- attr(M, "pointers")
  # hashes <- attr(M, "hashes")
  # kcounts <- attr(M, "kcounts") # NULL if not weighted

  # seqlengths <- apply(kcounts, 1, sum) + k - 1
  # minsql <- which.min(seqlengths)
  # scaledkc <- round(kcounts*(seqlengths[minsql]/seqlengths))
  kfreqs <- kcounts/apply(kcounts, 1, sum)

  # seqlengths <- attr(M, "seqlengths")
  # if(is.null(rownames(M))) rownames(M) <- paste0("S", 1:nrow(M))
  # labs <- rownames(M)
  # nuseq <- sum(!duplicates)
  # M <- M[!duplicates, ] # removes attrs
  #seqlengths <- sapply(x[!duplicates], length)
  ## initialize the tree
  tree <- 1
  attr(tree, "leaf") <- TRUE
  attr(tree, "sequences") <- 1:nuseq
  attr(tree, "height") <- 0
  # define recursive splitting functions
  topdown1 <- function(tree, d){ # d is the mbed matrix
    tree <- topdown2(tree, d = d)
    if(is.list(tree)) tree[] <- lapply(tree, topdown1, d = d)
    return(tree)
  }
  topdown2 <- function(node, d){
    if(!is.list(node) & length(attr(node, "sequences")) > 1){
      # fork leaves only
      seqs <- d[attr(node, "sequences"), , drop = FALSE]
      errfun <- function(er){# use when > 3 uniq hashes but kmeans throws error
        out <- list()
        nrs <- nrow(seqs)
        cls <- rep(1, nrs)
        cls[sample(1:nrs, 1)] <- 2 # peel one off
        out$cluster <- cls
        out$centers <- rbind(apply(seqs[cls == 1, , drop = FALSE], 2, mean),
                             apply(seqs[cls == 2, , drop = FALSE], 2, mean))
        return(out)
      }
      # membership <- tryCatch(kmeans(seqs, centers = 2)$cluster,
      #                 error = errfun, warning = errfun)
      km <- tryCatch(kmeans(seqs, centers = 2, ... = ...), error = errfun, warning = errfun)
      membership <- km$cluster
      centers <- km$centers
      tmpattr <- attributes(node)
      node <- vector(mode = "list", length = 2)
      attributes(node) <- tmpattr
      attr(node, "leaf") <- NULL
      attr(node, "avdist") <- sqrt(sum(abs(km$centers[1,] - km$centers[2,])^2))
      for(i in 1:2){
        node[[i]] <- 1
        attr(node[[i]], "height") <- attr(node, "height") - 1
        attr(node[[i]], "leaf") <- TRUE
        attr(node[[i]], "sequences") <- attr(node, "sequences")[membership == i]
        attr(node[[i]], "avdist") <- 0 # recalculated later
      }
    }
    return(node)
  }
  #  build tree recursively
  # tree <- topdown1(tree, d = scaledkc) # was M
  tree <- topdown1(tree, d = kfreqs) # was M
  tree <- remidpoint(tree)
  class(tree) <- "dendrogram"
  # model the branch lengths
  # if(weighted){
    # if(is.null(kcounts)) stop("Can't weight tree without kcounts")
    # avdist <- function(node, kcounts){
    #   if(is.list(node)){
    #     node1seqs <- attr(node[[1]], "sequences")
    #     node2seqs <- attr(node[[2]], "sequences")
    #     ln1s <- length(node1seqs)
    #     ln2s <- length(node2seqs)
    #     if(ln1s > 19) node1seqs <- sample(node1seqs, size = ceiling(log(ln1s, 2)^2))
    #     if(ln2s > 19) node1seqs <- sample(node2seqs, size = ceiling(log(ln2s, 2)^2))
    #       dists <- .kdist(kcounts, from = node1seqs - 1,
    #                       to = node2seqs - 1, seqlengths = seqlengths, k = k)
    #       attr(node, "avdist") <- mean(dists)
    #     # }
    #   }else{
    #     attr(node, "avdist") <- 0
    #   }
    #   return(node)
    # }
    # tree <- dendrapply(tree, avdist, kcounts = kcounts)
    reheight <- function(node){
      if(is.list(node)){
        node1edge <- max(0.0001, (attr(node, "avdist") - attr(node[[1]], "avdist"))/2)
        node2edge <- max(0.0001, (attr(node, "avdist") - attr(node[[2]], "avdist"))/2)
        attr(node[[1]], "height") <- attr(node, "height") - node1edge
        attr(node[[2]], "height") <- attr(node, "height") - node2edge
      }
      attr(node, "avdist") <- NULL # tidys up by removing attrs
      return(node)
    }
    reheight1 <- function(node){
      node <- reheight(node)
      if(is.list(node)) node[] <- lapply(node, reheight1)
      return(node)
    }
    tree <- reheight1(tree)
    class(tree) <- "dendrogram"
    tree <- reposition(tree)
  # }else{
  #   tree <- reposition(tree)
  #   tree <- ultrametricize(tree)
  # }
  if(any(duplicates)){
    reduplicate <- function(node, pointers){
      attr(node, "sequences") <- which(pointers %in% attr(node, "sequences"))
      if(is.leaf(node)){
        lams <- length(attr(node, "sequences"))
        if(lams > 1){
          labs <- attr(node, "label")
          hght <- attr(node, "height")
          seqs <- attr(node, "sequences")
          node <- vector(mode = "list", length = lams)
          attr(node, "height") <- hght
          attr(node, "sequences") <- seqs
          for(i in 1:lams){
            node[[i]] <- 1
            attr(node[[i]], "height") <- hght
            attr(node[[i]], "label") <- labs[i]
            attr(node[[i]], "sequences") <- seqs[i]
            attr(node[[i]], "leaf") <- TRUE
          }
        }
      }
      return(node)
    }
    tree <- dendrapply(tree, reduplicate, pointers)
    tree <- remidpoint(tree)
  }
  label <- function(node, labs){
    if(is.leaf(node)) attr(node, "label") <- labs[attr(node, "sequences")]
    return(node)
  }
  tree <- dendrapply(tree, label, labs = catchnames)
  rmseqs <- function(node){
    if(is.leaf(node)){
      tmpattr <- attributes(node)
      node[] <- tmpattr$sequences
      tmpattr$sequences <- NULL
      attributes(node) <- tmpattr
    }else{
      attr(node, "sequences") <- NULL
    }
    return(node)
  }
  tree <- dendrapply(tree, rmseqs)
  return(tree)
}
################################################################################

# # first temporarily remove duplicated sequences
# hashes <- sapply(x, function(s) paste(openssl::md5(as.vector(s))))
# duplicates <- duplicated(hashes)
# hasduplicates <- any(duplicates)
# if(hasduplicates){
#   pointers <- integer(length(x))
#   dhashes <- hashes[duplicates]
#   ndhashes <- hashes[!duplicates]
#   pointers[!duplicates] <- seq_along(ndhashes)
#   pointers[duplicates] <- sapply(dhashes, match, ndhashes)
#   fullx <- x
#   x <- x[!duplicates]
# }
# # then after reheighting tree:
# if(hasduplicates){
#   add_duplicates <- function(node, pointers){
#     attr(node, "nunique") <- length(attr(node, "sequences"))
#     attr(node, "sequences") <- which(pointers %in% attr(node, "sequences"))
#     attr(node, "ntotal") <- length(attr(node, "sequences"))
#     return(node)
#   }
#   #if(!quiet) cat("Repatriating duplicate sequences with tree\n")
#   tree <- dendrapply(tree, add_duplicates, pointers)
#   forknode2 <- function(node){
#     if(!is.list(node)){# fork leaves only
#       lns <- length(attr(node, "sequences"))
#       if(lns == 1) return(node)
#       membership <- c(1, rep(2, lns - 1))
#       tmpattr <- attributes(node)
#       node <- vector(mode = "list", length = 2)
#       attributes(node) <- tmpattr
#       attr(node, "leaf") <- NULL
#       for(i in 1:2){
#         node[[i]] <- 1
#         attr(node[[i]], "height") <- attr(node, "height") - 0.0001
#         attr(node[[i]], "leaf") <- TRUE
#         attr(node[[i]], "sequences") <- attr(node, "sequences")[membership == i]
#       }
#     }
#     return(node)
#   }
#   forknode1 <- function(node){
#     node <- forknode2(node)
#     if(is.list(node)) node[] <- lapply(node, forknode1)
#     return(node)
#   }
#   tree <- forknode1(tree)
#   tree <- remidpoint(tree)
#   x <- fullx
# }

#
# topdown2 <- function(node, d){
#   if(!is.list(node)){ # fork leaves only
#     seqs <- d[attr(node, "sequences"), , drop = FALSE]
#     nseq <- nrow(seqs)
#     allident <- FALSE
#     if(nseq == 1) return(node)
#     if(nseq == 2){
#       membership <- 1:2
#       if(identical(d[1,], d[2,])) allident <- TRUE
#     }else{
#       if(all(apply(seqs, 2, function(v) all(v[2:nseq] == v[1])))){
#         allident <- TRUE
#         membership <- c(1, rep(2, nseq - 1))
#       }else{
#         #test <<- seqs
#         membership <- kmeans(seqs, centers = 2)$cluster
#       }
#     }
#     tmpattr <- attributes(node)
#     node <- vector(mode = "list", length = 2)
#     attributes(node) <- tmpattr
#     attr(node, "leaf") <- NULL
#     for(i in 1:2){
#       node[[i]] <- 1
#       attr(node[[i]], "height") <- attr(node, "height") - if(allident) 0.0001 else 1
#       attr(node[[i]], "leaf") <- TRUE
#       attr(node[[i]], "sequences") <- attr(node, "sequences")[membership == i]
#     }
#   }
#   return(node)
# }
#
#
#
# #
# rdname topdown
# ################################################################################
# topdown.default <- function(x){
#   # x is a matrix with one row for each sequence
#   # each column is a dimension eg count for tuple k, dist to seed r, etc
#   # produces a dendogram via recursive k-means clustering
#   dupes <- duplicated.array(x, MARGIN = 1)
#   hasdupes <- any(dupes)
#   if(hasdupes){
#     whichdupes <- which(dupes)
#     lwd <- length(whichdupes)
#     tmp <- rbind(x[dupes, ], x)
#     dupes2 <- duplicated.array(tmp, MARGIN = 1)
#     whichdupes2 <- which(dupes2)
#     whichdupes2 <- whichdupes2[whichdupes2 > lwd] - lwd
#     m <- as.matrix(dist(x[whichdupes2, ]))
#     test <- as.dendrogram(hclust(dupdist, "average"))
#     test2 <- which(dupdist == 0)
#     zerodists <- which(m == 0 & lower.tri(m), arr.ind = T)
#   }
#   tree <- 1
#   attr(tree, "leaf") <- TRUE
#   attr(tree, "sequences") <- 1:nrow(x)
#   attr(tree, "height") <- 0
#   tree <- topdown1(tree, d = x)
#   tree <- remidpoint(tree)
#   class(tree) <- "dendrogram"
#   tree <- reposition(tree)
#   label <- function(node, d){
#     if(is.leaf(node)) attr(node, "label") <- rownames(d)[attr(node, "sequences")]
#     return(node)
#   }
#   tree <- dendrapply(tree, label, d = x)
#   return(tree)
# }


# topdown2 <- function(node, d){# d is the mbed object
#   if(!is.list(node) & length(attr(node, "sequences")) > 1){
#     # fork leaves only
#     seqs <- d[attr(node, "sequences"), , drop = FALSE]
#     hashes <- attr(d, "hashes")[attr(node, "sequences")]
#     nuniq <- 1
#     counter <- 2
#     uhashes <- character(2)
#     uhashes[1] <- hashes[1]
#     repeat{ # this is faster than length(unique(hashes)) for large nseq
#       if(hashes[counter] != uhashes[1] & hashes[counter] != uhashes[2]){
#         nuniq <- nuniq + 1
#         if(nuniq > 2) break
#         uhashes[2] <- hashes[counter]
#       }
#       if(counter == length(hashes)) break
#       counter <- counter + 1
#     }
#     errfun <- function(hashes){# use when > 3 uniq hashes but kmeans throws error
#       out <- rep(1, length(hashes))
#       out[hashes == sample(hashes, 1)] <- 2
#       return(out)
#     }
#     if(nuniq == 3){
#       membership <- tryCatch(kmeans(seqs, centers = 2)$cluster,
#                              error = errfun, warning = errfun)
#       nbranches <- 2
#       allid <- FALSE
#     }else if(nuniq == 2){
#       membership <- errfun(hashes)
#       nbranches <- 2
#       allid <- FALSE
#     }else{
#       nbranches <- nrow(seqs) # all identical
#       membership <- 1:nbranches
#       allid <- TRUE
#     }
#     tmpattr <- attributes(node)
#     node <- vector(mode = "list", length = nbranches)
#     attributes(node) <- tmpattr
#     attr(node, "leaf") <- NULL
#     for(i in 1:nbranches){
#       node[[i]] <- 1
#       attr(node[[i]], "height") <- attr(node, "height") - if(allid) 0 else 1
#       attr(node[[i]], "leaf") <- TRUE
#       attr(node[[i]], "sequences") <- attr(node, "sequences")[membership == i]
#     }
#   }
#   return(node)
# }
