#-------------------------------------------------------------------------------
# interpolate_dut1 
# 
# spline fit for Earth Orientation Parameters (EOP)
#
# data from http://maia.usno.navy.mil/
#
# @param utc datetime in utc
# @param fit_var variable to interpolate (ddt, x, y, dx, dy, lod, utc1_utc)
#
# @return seconds difference  32.184s + (TAI - UTC) - (UT1 - UTC)
# 
# @keywords internal
# 
#-------------------------------------------------------------------------------
interpolate_dut1 <- function(utc, fit_var, eop){
  
  # eop included in package, needs to be updated as new leap seconds are added
  #approx(eop[['datetime']], eop[[fit_var]], utc, ties = 'ordered')$y
  stats::spline(eop[['datetime']], eop[[fit_var]],
                xout = utc, ties = 'ordered')$y
  
  
}

