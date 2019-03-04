#-------------------------------------------------------------------------------
#' delta_t
#'
#' delta_t returns the delta t
#' data from http://maia.usno.navy.mil/
#'
#' @param utc datetime in utc
#'
#' @return seconds difference  32.184s + (TAI - UTC) - (UT1 - UTC)
#'
#-------------------------------------------------------------------------------
delta_t <- function(utc){
  
  # dut1 included in package, needs to be updated as new leap seconds are added
  
  #approx(dut1$datetime, dut1$ddt, utc, ties = 'ordered')$y
  stats::spline(dut1$datetime, dut1$ddt, xout = utc, ties = 'ordered')$y
  
  
}

