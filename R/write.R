#' Write a dendrogram object to parenthetic text.
#'
#' \code{write.dendrogram} writes a dendrogram object to a text file
#'   in Newick (New Hampshire) format.
#'
#' @param x an object of class \code{"dendrogram"}.
#' @param file a character string naming a file or connection to write the
#'   output to. If no file path is specified or \code{file = ""} the result
#'   will be printed to the console.
#' @param append a logical value indicating whether the output should be
#'   appended to the file. If \code{append = FALSE} the contents of the
#'   file will be overwritten (the default setting).
#' @param edges a logical value indicating whether edge weights should be
#'   included in the output string.
#' @param ... further arguments to be passed to \code{format}.Used to
#'   specify the numbering style of the edge weights (if edges = TRUE).
#' @seealso \code{\link{read.dendrogram}} to create a \code{"dendrogram"}
#'   object from a text file.
#'   The \code{\link[ape]{write.tree}} function in the \code{\link[ape]{ape}}
#'   package performs a similar operation for objects of class \code{"phylo"}
#'   and \code{"multiPhylo"}.
#' @examples
#'   arrests.hclust <- hclust(dist(USArrests[1:6,]), "average")
#'   arrests.dendro <- as.dendrogram(arrests.hclust)
#'   write.dendrogram(arrests.dendro, digits = 3)
#'
################################################################################
write.dendrogram <- function(x, file = "", append = FALSE, edges = TRUE, ...){
  if(!(inherits(x, "dendrogram"))) stop("Input object must be of class
                                        'dendrogram'")
  renameLeaves <- function(y){
    if(is.leaf(y)){
      tmp <- attr(y, "label")
      if(grepl("[^abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789]",
               tmp)){
        tmp <- paste0(c("'", tmp, "'"), collapse = "")
      }
      y[1] <- tmp
    }
    y
  }
  x <- dendrapply(x, renameLeaves)
  xnames <- unlist(x)
  ynames <- paste0("S", seq_along(xnames) + 10000)
  renameLeaves2 <- function(y){
    if(is.leaf(y)){
      y[1] <- ynames[match(y[1], xnames)]
    }
    y
  }
  x <- dendrapply(x, renameLeaves2)
  x <- dendrapply(x, unclass)
  addEdges <- function(y){
    if(is.list(y)){
      y[] <- lapply(y, function(z){
        attr(z, "edge") <- format(attr(y, "height") - attr(z, "height"),
                                  ... = ...)#scientific = FALSE)
        z
      })
      attributes(y)[names(attributes(y)) != "edge"] <- NULL
      y[] <- lapply(y, addEdges)
    }else{
      attributes(y)[names(attributes(y)) != "edge"] <- NULL
    }
    y
  }
  x <- addEdges(x)
  attr(x, "edge") <- 0
  tmp <- deparse(x)
  tmp <- paste0(tmp, collapse = "")
  tmp <- gsub(" ", "", tmp)
  tmp <- gsub("edge=\"([-0-9Ee.]+)\"", "edge=\\1", tmp)
  #for(i in seq_along(xnames)) tmp <- gsub(xnames[i], ynames[i], tmp)
  tmp2 <- gsub("list", ";", tmp) # needs to be a special single char - alternative?
  while(grepl("structure", tmp2)){
    tmp2 <- gsub("structure\\(\"([^\"]*)\",edge=([-0-9Ee.]+)\\)", "\\1:\\2", tmp2)
    tmp2 <- gsub(";\\(([^;\\)]*)\\)", "\"openbracket\\1closebracket\"", tmp2)
  }
  tmp3 <- gsub("openbracket", "\\(", tmp2)
  tmp3 <- gsub("closebracket", "\\)", tmp3)
  res <- gsub("(.*):0$", "\\1;", tmp3)
  if(!edges){
    res <- gsub(":[-0-9Ee.]+", "", res)
  }
  for(i in seq_along(xnames)) res <- gsub(ynames[i], xnames[i], res)
  if(file == ""){
    return(res)
  }else{
    cat(res, file = file, append = append, sep = "\n")
  }
}



