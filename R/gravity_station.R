#-------------------------------------------------------------------------------
# gravity_station
#
# \code{gravity_station} returns the estimated gravity for a latitude.
# http://the-mostly.ru/misc/local_gravity_online_calculator.html
#
# @param latitude the station latitude 
# @param elevation the station elevation in meters
#
# @return gravity for location in m/s^2
# 
# @keywords internal
# 
#-------------------------------------------------------------------------------
gravity_station <- function(latitude, elevation){
  
  to_radians <- pi/180
  
  a <- 6371000
  sin_lat <- sin(latitude * to_radians)^2

  9.7803267714 * (1 + 0.00193185138639 * sin_lat) /
    sqrt(1 - 0.00669437999013 * sin_lat) * 
    (1 + elevation / a)^(-2) 
  
  # eccen  <- 6.69439795140e-3
  # 
  # 9.78032677 * (1.0 + 0.001931851353 * sin_lat^2) /
  #   sqrt(1.0 - eccen * sin_lat^2) - 0.3086e-5 * elevation

}