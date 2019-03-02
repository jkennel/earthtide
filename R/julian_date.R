#' utc_julian
#'
#' \code{utc_julian} returns the julian datetime from UTC
#'
#' @param utc datetime in utc
#'
#' @return julian datetime
#'
utc_julian <- function(utc) {
  
  (as.numeric(utc) / 86400.0) + 2440587.5;
  
}


#' utc_julian_2000
#'
#' \code{utc_julian_2000} returns the j2000 datetime from UTC
#'
#' @param utc datetime in utc
#'
#' @return julian datetime
#'
utc_julian_2000 <- function(utc) {
  
  (utc_julian(utc) - 2451545.0) / 36525.0

}

#' utc_mod_julian
#'
#' \code{utc_mod_julian} returns the modified julian datetime from UTC
#'
#' @param utc datetime in utc
#'
#' @return julian datetime
#'
utc_mod_julian <- function(utc) {
  
  (utc_julian(utc) - 2400000.5)
  
}


