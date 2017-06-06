#' Count k-mers.
#'
#' Count all k-letter words in a sequence or set of sequences by
#'   sliding a window of length k over the sequence(s).
#'
#' @param x a matrix of aligned sequences, a list of unaligned sequences,
#'   or a vector representing a single sequence.
#'   Accepted modes are "character" and "raw" (the latter is
#'   for "DNAbin" and "AAbin" objects).
#' @param k integer representing the k-mer size. Defaults to 5.
#'   Note that high values of k
#'   may be slow to compute and use a lot of memory due to the large numbers
#'   of calculations required, particularly when the residue alphabet is
#'   also large.
#' @param residues either NULL (default; the residue alphabet is automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Specifying the residue alphabet is therefore
#'   recommended unless x is a "DNAbin" or "AAbin" object.
#' @param gap the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} or \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @return a matrix of k-mer counts with one row for each sequence and n^k
#'   columns (where n is the size of the residue alphabet and k is the
#'   k-mer size)
#' @details TBA
#' @author Shaun Wilkinson
#' @references TBA
#' @seealso \code{\link{kdistance}} for distance matrices based on k-mer counting
#' @examples
#'   ## compute a matrix of k-mer counts for the woodmouse data (ape package)
#'   ## with a k-mer size of 3
#'   library(ape)
#'   data(woodmouse)
#'   x <- kcount(woodmouse, k = 3)
#'   ## 15 row matrix with 64 columns for nucleotide 3-mers AAA, AAC, ... TTT
#'   ##
#'   ## convert to AAbin object and repeat the operation
#'   y <- kcount(ape::trans(woodmouse, 2), k = 2)
#'   ## 15 row matrix with 400 columns for amino acid 2-mers AA, AB, ... , YY
################################################################################
kcount <- function(x, k = 5, residues = NULL, gap = "-"){
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  residues <- .alphadetect(x, residues = residues, gap = gap)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  if(is.matrix(x)) x <- .unalign(x, gap = gap)
  if(!is.list(x)) x <- list(x)
  nseq <- length(x)
  if(is.null(names(x))) names(x) <- paste0("S", 1:nseq)
  seqalongx <- seq_along(x)
  if(DNA){
    arity <- 4
    x <- lapply(x, function(s) s[!(s %in% as.raw(c(2, 4, 240)))])
    seqlengths <- sapply(x, length)
    if(min(seqlengths) < k) stop("minimum sequence length is less than k")
    kcounts <- .kcountDNA(x, k = k)
  }else{
    tuplecount <- function(y, k, arity){
      tuplemat <- matrix(nrow = k, ncol = length(y) - k + 1)
      for(i in 1:k) tuplemat[i, ] <- y[i:(length(y) - (k - i))]
      res <- apply(tuplemat, 2, .decimal, from = arity) + 1
      res <- tabulate(res, nbins = arity^k)
      return(res)
    }
    if(AA){
      arity <- if(k > 3) 6 else 20 # compress AA alphabet for high k values
      x <- .encodeAA(x, arity = arity, na.rm = TRUE)
    }else{
      #residues <- .alphadetect(x, residues = residues, gap = gap)
      arity <- length(residues)
      if(k > 2 & arity >= 20) stop("Unable to calculate distance matrix for
                               large k and large alphabet size. If residues
                               are amino acids consider converting to AAbin
                               object for compression")
      modes <- lapply(x, mode)
      if(!(all(modes == "integer") | all(modes == "numeric"))){
        x <- .encodeCH(x, residues = residues, na.rm = TRUE)
      }
    }
    seqlengths <- sapply(x, length)
    if(min(seqlengths) < k) stop("minimum sequence length is less than k")
    kcounts <- t(sapply(x, tuplecount, k, arity))
  }
  # label columns with k-mer words
  indices <- matrix(1L, nrow = k, ncol = arity^k)
  ntimes <- 1
  counter <- k - 1
  for(i in 1:k){
    indices[i, ] <- rep(rep(1:arity, each = arity^counter), times = ntimes)
    ntimes <- ntimes * arity
    counter <- counter - 1
  }
  colnames(kcounts) <- apply(indices, 2, function(i) paste0(residues[i], collapse = ""))
  return(kcounts)
}
################################################################################
#' K-mer distance matrix computation.
#'
#' Computes the matrix of k-mer distances between all pairwise comparisons
#'   of a set of sequences.
#'
#' @param x a matrix of aligned sequences or a list of unaligned sequences.
#'   Accepted modes are "character" and "raw" (the latter is
#'   for "DNAbin" and "AAbin" objects).
#' @param k integer representing the k-mer size to be used for calculating
#'   the distance matrix. Defaults to 5. Note that high values of k
#'   may be slow to compute and use a lot of memory due to the large numbers
#'   of calculations required, particularly when the residue alphabet is
#'   also large.
#' @param method a character string giving the k-mer distance measure
#'   to be used. Currently the available options are \code{"edgar"} (default;
#'   see Edgar (2004) for details) and the standard methods avaialable for
#'   the base function "dist" ("euclidean", "maximum", "manhattan", "canberra",
#'   "binary" and "minkowski").
#' @param residues either NULL (default; the residue alphabet is automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Specifying the residue alphabet is therefore
#'   recommended unless x is a "DNAbin" or "AAbin" object.
#' @param gap the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} or \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @param ... further arguments to be passed to \code{"as.dist"}.
#' @return a distance matrix of class \code{"dist"}
#' @details this function computes the N x N k-mer distance matrix,
#'   where N is the number of sequences in the dataset.
#' @author Shaun Wilkinson
#' @references
#'   Edgar RC (2004) Local homology recognition and distance measures in
#'   linear time using compressed amino acid alphabets.
#'   \emph{Nucleic Acids Research}, \strong{32}, 380-385.
#'
#' @seealso \code{\link{mbed}} for leaner distance matrices
#' @examples
#'   ## compute the k-mer distance for the woodmouse dataset in the ape package
#'   ## with a k-mer size of 5
#'   library(ape)
#'   data(woodmouse)
#'   ## trim gappy ends for global alignment
#'   woodmouse <- woodmouse[, apply(woodmouse, 2, function(v) !any(v == 0xf0))]
#'   woodmouse.dist <- kdistance(woodmouse, k = 5)
#'   ## cluster and plot UPGMA tree
#'   woodmouse.tree <- as.dendrogram(hclust(woodmouse.dist, "average"))
#'   plot(woodmouse.tree)
################################################################################
kdistance <- function(x, k = 5, method = "edgar", residues = NULL,
                      gap = "-", ...){
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  residues <- .alphadetect(x, residues = residues, gap = gap)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  if(is.matrix(x)) x <- .unalign(x, gap = gap)
  nseq <- length(x)
  if(is.null(names(x))) names(x) <- paste0("S", 1:nseq)
  seqalongx <- seq_along(x)
  kcounts <- kcount(x, k = k, residues = residues, gap = gap)
  seqlengths <- apply(kcounts, 1, sum) + k - 1
  ## not sapply(x, length) in case of gaps, unknowns etc
  ## which will be picked up by kcount
  if(method == "edgar"){
    d <- .kdist(kcounts, from = seqalongx - 1, to = seqalongx - 1,
                        seqlengths = seqlengths, k = k)
  }else{
    freqs <- kcounts/(seqlengths - k + 1)
    d <- dist(freqs, method = method)
  }
  return(as.dist(d, ... = ...))
}
################################################################################
#' Convert sequences to vectors of distances to \emph{n} seed sequences.
#'
#' The \code{mbed} function takes a list of sequences and returns a matrix of
#'   distances to a subset of seed sequences using the method outlined
#'   in Blacksheilds et al. (2010).
#'
#' @param x a matrix of aligned sequences or a list of unaligned sequences.
#'   Accepted modes are "character" and "raw" (the latter is for "DNAbin"
#'   and "AAbin" objects).
#' @param k integer representing the k-mer size to be used for calculating
#'   the distance matrix. Defaults to 5. Note that high values of k
#'   may be slow to compute and use a lot of memory due to the large numbers
#'   of calculations required, particularly when the residue alphabet is
#'   also large.
#' @param seeds optional integer vector indicating which sequences should
#'   be used as the seed sequences. If \code{seeds = NULL} a set of
#'   log(N, 2)^2 non-identical sequences is randomly selected from the
#'   sequence set (where N is the number of sequences; see Blacksheilds et al.
#'   2010). Alternatively, if \code{seeds = 'all'} a standard N x N
#'   distance matrix is computed.
#' @param residues either NULL (default; emitted residues are automatically
#'   detected from the sequences), a case sensitive character vector
#'   specifying the residue alphabet, or one of the character strings
#'   "RNA", "DNA", "AA", "AMINO". Note that the default option can be slow for
#'   large lists of character vectors. Specifying the residue alphabet is therefore
#'   recommended unless x is a "DNAbin" or "AAbin" object.
#' @param gap the character used to represent gaps in the alignment matrix
#'   (if applicable). Ignored for \code{"DNAbin"} or \code{"AAbin"} objects.
#'   Defaults to "-" otherwise.
#' @param counts logical indicating whether the (usually large) matrix of
#'   k-mer counts should be returned as an attribute of the returned
#'   object. Defaults to FALSE to avoid excessive memory usage.
#' @return returns a N x log(N, 2)^2 matrix of class "mbed" (where N is the
#'   number of sequences). The returned
#'   object has additional attributes including a the 'seeds' vector and
#'   a matrix of k-mer counts with one row for each sequence and 4^k columns.
#' @details This function deals with ambiguities by assigning counts proportionally.
#'   For example the motif ACRTG would assign the 5-mers ACATG and ACGTG counts of
#'   0.5 each. This algorithm is O N * 4^k in memory and time complexity so can be very
#'   slow and memory hungry for larger values of k. Further details of the embedding
#'   process can be found in Blacksheilds et al. (2010). The k-mer distance measure
#'   used is that of Edgar (2004).
#' @author Shaun Wilkinson
#' @references
#'   Blackshields G, Sievers F, Shi W, Wilm A, Higgins DG (2010) Sequence embedding
#'   for fast construction of guide trees for multiple sequence alignment.
#'   \emph{Algorithms for Molecular Biology}, \strong{5}, 21.
#'   Edgar RC (2004) Local homology recognition and distance measures in
#'   linear time using compressed amino acid alphabets.
#'   \emph{Nucleic Acids Research}, \strong{32}, 380-385.
#' @seealso \code{\link{kdistance}} for full N x N distance matrix computation
#' @examples
#'   library(ape)
#'   data(woodmouse)
#'   ## randomly select three sequences as seeds
#'   set.seed(999)
#'   seeds <- sample(1:15, size = 3)
#'   ## embed the woodmouse dataset in three dimensions
#'   woodmouse.mbed <- mbed(woodmouse, seeds = seeds, k = 5)
#'   ## print the distance matrix without attributes
#'   print(woodmouse.mbed[,], digits = 2)
################################################################################
mbed <- function(x, seeds = NULL, k = 5, residues = NULL, gap = "-",
                 counts = FALSE){
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
  residues <- .alphadetect(x, residues = residues, gap = gap)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  if(is.matrix(x)) x <- .unalign(x, gap = gap)
  nseq <- length(x)
  if(is.null(names(x))) names(x) <- paste0("S", 1:nseq)
  if(!is.null(seeds)){
    if(identical(seeds, "all")) seeds <- seq_along(x)
    stopifnot(mode(seeds) %in% c("numeric", "integer"),
              max(seeds) <= nseq,
              min(seeds) > 0)
  }
  hashes <- .digest(x, simplify = TRUE)
  # hashes <- sapply(x, function(s) paste(openssl::md5(as.vector(s))))
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
    catchnames <- names(x)
    x <- x[!duplicates]
  }else{
    pointers <- seq_along(x)
  }
  seqalongx <- seq_along(x)
  if(DNA){
    x <- lapply(x, function(s) s[!(s %in% as.raw(c(2, 4, 240)))])
    seqlengths <- sapply(x, length)
    if(min(seqlengths) < k) stop("minimum sequence length is less than k")
    kcounts <- .kcountDNA(x, k = k)
  }else{
    tuplecount <- function(y, k, arity){
      tuplemat <- matrix(nrow = k, ncol = length(y) - k + 1)
      for(i in 1:k) tuplemat[i, ] <- y[i:(length(y) - (k - i))]
      res <- apply(tuplemat, 2, .decimal, from = arity) + 1
      res <- tabulate(res, nbins = arity^k)
      return(res)
    }
    if(AA){
      arity <- if(k > 2) 6 else 20 # compress AA alphabet for high k values
      x <- .encodeAA(x, arity = arity, na.rm = TRUE)
    }else{
      arity <- length(residues)
      if(k > 2 & arity >= 20) stop("Unable to calculate distance matrix for
                               large k and large alphabet size. If residues
                               are amino acids consider converting to AAbin
                               object for automatic compression, else reduce k")
      modes <- lapply(x, mode)
      if(!(all(modes == "integer") | all(modes == "numeric"))){
        x <- .encodeCH(x, residues = residues, na.rm = TRUE)
      }
    }
    seqlengths <- sapply(x, length)
    if(min(seqlengths) < k) stop("minimum sequence length is less than k")
    kcounts <- t(sapply(x, tuplecount, k, arity))
  }
  if(is.null(seeds)){
    nseeds <- ceiling(log(nuseq, 2)^2)
    allseeds <- nseeds >= nuseq
    if(allseeds){
      nseeds <- nuseq
      seeds <- seqalongx
    }else{
      suppressWarnings(groups <- kmeans(kcounts, centers = nseeds)$cluster)
      seeds <- match(1:nseeds, groups)
    }
    ## LLR algorithm see Blacksheilds et al. 2010
    # seeds <- sort(sample(seqalongx, size = nseeds))
  }else{
    seeds <- unique(pointers[seeds])
    nseeds <- length(seeds)
  }
  res <- .kdist(kcounts, from = seqalongx - 1, to = seeds - 1,
                seqlengths = seqlengths, k = k)
  if(any(duplicates)){
    tmp <- matrix(nrow = nseq, ncol = ncol(res))
    rownames(tmp) <- catchnames
    colnames(tmp) <- names(x)[seeds]
    tmp[!duplicates, ] <- res
    #for(i in which(duplicates)) tmp[i, ] <- res[pointers[i], ]
    tmp[duplicates, ] <- res[pointers[duplicates], ]
    res <- tmp
    # also refll kcounts
    if(counts){
      tmpkc <- matrix(nrow = nseq, ncol = ncol(kcounts))
      rownames(tmpkc) <- catchnames
      tmpkc[!duplicates, ] <- kcounts
      tmpkc[duplicates, ] <- kcounts[pointers[duplicates], ]
      kcounts <- tmpkc
    }
  }
  attr(res, "seeds") <- seeds # integer vector
  if(counts) attr(res, "kcounts") <- kcounts
  rm(kcounts)
  gc()
  attr(res, "duplicates") <- duplicates
  attr(res, "pointers") <- pointers
  attr(res, "hashes") <- hashes
  attr(res, "seqlengths") <- seqlengths
  class(res) <- "mbed"
  return(res)
}
################################################################################
