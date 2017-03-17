#' A package for viewing, editing and publishing phylogenetic trees in R.
#'
#' The phylogram package contains functions for reading, viewing and
#'   editing phylogenetic trees as deeply nested lists using R's
#'   'dendrogram' object type. This enables users to perform both
#'   top-down and bottom-up recursive tree operations such as splitting
#'   and neighbor joining. For compatibility with other programs and
#'   packages, trees can be imported and exported
#'   in the 'Newick' or 'New Hampshire' text format.
#'
#' @section phylogram functions:
#'   The \code{read.phylogram} and \code{write.phylogram} functions import
#'   and export phylogram objects as newick text strings via a file or
#'   connection.
#'   \code{prune} and \code{ladder} are tree editing functions, that are
#'   used to remove branches (based on regular expression pattern matching)
#'   and reorder branches, respectively. Other editing functions include
#'   \code{remidpoint}, \code{reposition}, and \code{ultrametricize}.
#'
#' @docType package
#' @name phylogram
################################################################################
NULL
