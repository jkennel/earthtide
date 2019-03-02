#' ut1_utc
#'
#' \code{tai_ut1} gets seconds difference UT1-UTC (valid after 1992)
#' 
#' @param utc datetime in utc
#'
#' @return seconds difference UT1-UTC
#'
ut1_utc <- function(utc){
  
  approx(dut1$datetime, dut1$ut1_utc, utc, ties = 'ordered')$y
  
}


