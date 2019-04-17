#' earthtide: R port of the earth tide processing package ETERNA 
#' (by Hans-Georg Wenzel) including the Kudryavtsev wave catalog. 
#'
#' The goal of this package is to generate synthetic earth tides for use in the
#' R programming language and in particular environmental models. Code was 
#' parallized and refactored to minimize duplication, and to allow for future 
#' improvements.
#'
#'
#' You can learn about the earthtide package in the vignettes:
#' \code{browseVignettes(package = "earthtide")}
#'
#'
#' @references Hartmann, T., Wenzel, H.-G., 1995. The HW95 tidal potential catalogue. Geophys. Res. Lett. 22, 3553-3556. \doi{10.1029/95GL03324}
#' @references Kudryavtsev, S.M., 2004. Improved harmonic development of the Earth tide-generating potential. J. Geod. 77, 829-838. \doi{10.1007/s00190-003-0361-2}
#' @references Wenzel, H.G., 1996. The nanogal software: Earth tide data processing package ETERNA 3.30. Bull. Inf. Marées Terrestres, 124, pp.9425-9439.  \url{http://www.eas.slu.edu/GGP/ETERNA34/MANUAL/ETERNA33.HTM}
#' 
#' @useDynLib earthtide, .registration = TRUE
#' @docType package
#' @aliases earthtide-package
#' @importFrom R6 R6Class
#' @importFrom stats approx
#' @importFrom utils data download.file
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
"_PACKAGE"
