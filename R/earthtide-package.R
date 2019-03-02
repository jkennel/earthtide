#' earthtide: R port of the earth tide processing package ETERNA 
#' (by Hans-Georg Wenzel) including the Kudryavtsev wave catalog. 
#'
#' The goal of this package is to generate synthetic earth tides for use in the
#' R programming language and in particular environmental models. Code was 
#' parallized and refactored to minimize duplication, and to allow for future 
#' improvements
#'
#'
#' You can learn about the earthtide package in the vignettes:
#' `browseVignettes(package = "earthtide")`
#'
#'
#' @useDynLib earthtide, .registration = TRUE
#' @docType package
#' @aliases earthtide-package
#' @importFrom R6 R6Class
#' @importFrom stats approx
#' @importFrom utils data download.file
#' @importFrom utils data download.file
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
"_PACKAGE"
