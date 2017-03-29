#' Print summary methods.
#' @param x object of various classes.
#' @param ... additional arguments to be passed between methods.
#' @return NULL (invisibly)
#' @author Shaun Wilkinson
################################################################################
print.mbed <- function(x, ...){
  cat("Matrix of", nrow(x),
      "embedded sequences, represented as vectors of distances to",
      ncol(x),
      "seed sequences")
}
################################################################################
