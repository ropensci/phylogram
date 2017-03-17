#' Read a dendrogram from parenthetic text.
#'
#' \code{read.dendrogram} parses a text file or character string in Newick
#'   (New Hampshire) format and creates an object of class \code{"dendrogram"}.
#'
#' @param file character string giving a valid path to the file from
#'   where to read the data.
#' @param text optional character string in lieu of a "file" argument.
#'   If a text argument is provided instead of a file path, the data
#'   are read via a text connection.
#' @param edges logical indicating whether edge weights
#'   provided in the Newick string should be retained in the returned
#'   object (defaults to TRUE).
#' @param ... further arguments to be passed to \code{scan}.
#' @details There are many varying interpretations of the Newick/
#'   New Hampshire file format. This function attempts to adhere to
#'   that outlined by Joe Felsenstein
#'   \href{http://evolution.genetics.washington.edu/phylip/newicktree.html}{here}.
#'   Newick strings with labeled inner nodes (for example
#'   "(B:6.0,(A:5.0,C:3.0,E:4.0)Ancestor1:5.0,D:11.0);")
#'   are accepted, however the current version removes the inner node
#'   labels. Comments enclosed in square brackets are also discarded.
#' @return an object of class \code{"dendrogram"}.
#' @author Shaun Wilkinson
#' @references
#'   \url{http://evolution.genetics.washington.edu/phylip/newicktree.html}
#'   \url{http://evolution.genetics.washington.edu/phylip/newick_doc.html}
#' @seealso
#'   \code{\link{write.dendrogram}} writes an object of
#'   class \code{"dendrogram"} to a text string.
#'   The \code{\link[ape]{read.tree}} function in the
#'   \code{\link[ape]{ape}} package performs a similar operation for objects
#'   of class \code{"phylo"} and \code{"multiPhylo"}.
#' @examples
#'   newick <- "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
#'   dendro <- read.dendrogram(text = newick)
#'   plot(dendro, horiz = TRUE)
#'
################################################################################
read.dendrogram <- function(file = "", text = NULL, edges = TRUE, ...){
  if(!is.null(text)){
    if(!is.character(text))
      stop("Argument 'text' must be of mode character.")
    x <- text
  }else{
    x <- scan(file = file, what = "", sep = "\n", quiet = TRUE, ...)
  }
  if (identical(x, character(0))) {
    warning("Empty character string.")
    return(NULL)
  }
  # collapse vector to a single string if necessary
  x <- paste0(x, collapse = "")
  # enclose entire string in brackets (excluding ;)
  xsplit <- strsplit(x, split = "")[[1]]
  if(xsplit[length(xsplit) - 1] != ")"){
    xsplit <- c("(", xsplit[-length(xsplit)], ");")
    x <- paste0(xsplit, collapse = "")
  }
  has.comments <- grepl("\\[", x) | grepl("\\]", x)
  if(has.comments){
    opens <- which(xsplit == "[")
    closes <- which(xsplit == "]")
    if(length(opens) != length(closes)) stop("Invalid metacharacters in Newick string")
    comments <- unlist(mapply(":", opens, closes))
    xsplit <- xsplit[-comments]
    x <- paste0(xsplit, collapse = "")
  }
  if(!edges){
    x <- gsub(":([-0-9Ee.]+)", "", x)
    has.edges <- FALSE
  }else{
    has.edges <- grepl(":", x)
  }
  has.unmatched.singlequotes <- grepl("''", x)
  if(has.unmatched.singlequotes) x <- gsub("''", "singlequote", x)
  # rectified later with fixnames function
  fun1 <- function(s, has.edges){ # a string, applied to odds
    # (not enclosed in single quotes)
    # Underscore characters outside unquoted labels are converted to blanks.
    #s <- gsub("_", "", s)
    s <- gsub(" ", "", s)
    # blank leaves are renamed, rectified later with fixnames function
    s <- gsub("\\( *,", "\\(unnamedleaf,", s)
    while(grepl(", *,", s)) s <- gsub(", *,", ",unnamedleaf,", s)
    s <- gsub(", *\\)", ",unnamedleaf\\)", s)
    # remove inner nodes (for now)
    s <- gsub("\\)[^,:;\\(\\)]+([,:;\\(\\)])", "\\)\\1", s)
    if(has.edges){
      s <- gsub("([\\(,])([^\\(\\),]+):", "\\1'\\2':", s)
      s <- gsub(";", ":1", s)
    }else{
      s <- gsub("([\\(,])([^\\(\\),]+)", "\\1\'\\2':1", s)
      s <- gsub("\\)","\\):1", s)
      s <- gsub(";", "", s)
    }
    return(s)
  }
  has.singlequotes <- grepl("'", x)
  if(has.singlequotes){
    tmp <- strsplit(x, split = "'")[[1]]
    evens <- seq(from = 2, to = length(tmp), by = 2)
    odds <- seq(from = 1, to = length(tmp), by = 2) # not names in single quotes
    tmp[odds] <- unname(sapply(tmp[odds], fun1, has.edges = has.edges))
    tmp[evens] <- unname(sapply(tmp[evens], function(s) paste0("'", s, "'", if(has.edges) "" else ":1")))
    tmp <- paste0(tmp, collapse = "")
  }else{
    tmp <- fun1(x, has.edges = has.edges)
  }
  tmp2 <- strsplit(tmp, split = "'")[[1]]
  evens <- seq(from = 2, to = length(tmp2), by = 2)
  leafnames <- tmp2[evens]
  has.unnamed.leaves <- any(grepl("unnamedleaf", leafnames))
  newleafnames <- paste0("L", seq(100001, 100000 + length(leafnames)))
  tmp2[evens] <- newleafnames
  tmp2 <- paste0(tmp2, collapse = "")
  tmp3 <- tmp2
  tree <- lapply(leafnames, function(e) e)
  names(tree) <- newleafnames
  innernodecount <- 100001
  while(grepl("\\([-LI0123456789Ee,:.]+\\)", tmp3)){
    # I is inner node, L is leaf, max 900000 leaves
    tojoin <- gsub(".*\\(([-LI0123456789Ee,:.]+)\\).*", "\\1", tmp3)
    tojoin2 <- strsplit(tojoin, split = ",")[[1]]
    ntojoin2 <- length(tojoin2)
    whichntj <- integer(ntojoin2)
    newnode <- vector(mode = "list", length = ntojoin2)
    for(s in seq_along(tojoin2)){
      nameedge <- strsplit(tojoin2[s], split = ":")[[1]]
      whichntj[s] <- match(nameedge[1], names(tree))
      attr(tree[[whichntj[s]]], "edge") <- as.numeric(nameedge[2])
      newnode[[s]] <- tree[[whichntj[s]]]
    }
    newnodename <- paste0("I", innernodecount)
    tree[[newnodename]] <- newnode
    tree <- tree[-whichntj]
    tmp3 <- gsub(paste0("\\(", tojoin, "\\)"), newnodename, tmp3)
    innernodecount <- innernodecount + 1
  }
  if(!grepl("^[LI][0-9]{6}:([-0-9Ee.]+)$", tmp3)) warning("Incomplete tree parse")
  #if(length(tree) == 1 & is.list(tree[[1]])) tree <- tree[[1]]
  tree <- tree[[1]]
  attr(tree, "edge") <- 0
  # convert nested list to dendrogram object by setting attributes recursvely
  setnodeattr <- function(x, leafnames){ # x is a nested list with 'edge' attributes, leafnmes is a character vector
    if(is.list(x)){
      cladesizes <- sapply(x, function(y) length(unlist(y)))
      attr(x, "members") <- sum(cladesizes)
      attr(x, "midpoint") <- if(length(cladesizes) > 1){
        ((cladesizes[1] - 1)/2 + (cladesizes[1] + (cladesizes[2] - 1)/2))/2
      }else 0
      if(is.null(attr(x, "height"))) attr(x, "height") <- 0
      setleafattr <- function(y, leafnames){
        attr(y, "height") <- attr(x, "height") - attr(y, "edge")
        attr(y, "edge") <- NULL
        if(!(is.list(y))){
          attr(y, "label") <- y[1]
          attr(y, "leaf") <- TRUE
          attr(y, "members") <- 1
          y[] <- match(y, leafnames)
          mode(y) <- "integer"
        }
        y
      }
      x[] <- lapply(x, setleafattr, leafnames = leafnames)
      x[] <- lapply(x, setnodeattr, leafnames = leafnames)
    }
    x
  }
  res <- setnodeattr(tree, leafnames = leafnames)
  if(!is.list(res)){
    attr(res, "height") <- attr(res, "edge")
    attr(res, "label") <- res[1]
    attr(res, "members") <- 1
    attr(res, "leaf") <- TRUE
  }
  attr(res, "edge") <- NULL
  # now convert to dendrogram
  attr(res, "class") <- "dendrogram"
  res <- reposition(res)
  # min.height <- min(unlist(dendrapply(res, attr, "height")))
  # reposition <- function(y, min.height){ # y is a dendrogram
  #   attr(y, "height") <- attr(y, "height") - min.height
  #   y
  # }
  # res <- dendrapply(res, reposition, min.height = min.height)
  fixnames <- function(y){
    if(!(is.list(y))){
      attr(y, "label") <- gsub("unnamedleaf", "", attr(y, "label"))
      attr(y, "label") <- gsub("singlequote", "'", attr(y, "label"))
    }
    y
  }
  if(has.unnamed.leaves | has.unmatched.singlequotes){
    res <- dendrapply(res, fixnames)
  }
  if(!has.edges) res <- ultrametricize(res)
  return(res)
}
################################################################################
