.prepare_station <- function(self, latitude, longitude, elevation,
                             azimuth, gravity, earth_radius,
                             earth_eccen) {
  
  angular_velocity = 7.292115e-5
  
  to_radians <- pi/180
  
  if (gravity == 0) {
    gravity <- gravity_station(latitude, elevation)
  } else {
    gravity <- gravity
  }
  
  sin_lat <- sin(latitude * to_radians)
  cos_lat <- cos(latitude * to_radians)
  sin_lon <- sin(longitude * to_radians)
  cos_lon <- cos(longitude * to_radians)
  

  # pole tide
  x_pol  <- interpolate_dut1(self$datetime$utc, 'x', self$datetime$eop)
  y_pol  <- interpolate_dut1(self$datetime$utc, 'y', self$datetime$eop) 
  dx_pol <- interpolate_dut1(self$datetime$utc, 'dx', self$datetime$eop) 
  dy_pol <- interpolate_dut1(self$datetime$utc, 'dy', self$datetime$eop) 
  
  self$pole_t <- 1.16 * 2.0 * angular_velocity^2 * 
    earth_radius * cos_lat * sin_lat * ((x_pol) * cos_lon - 
                                          (y_pol) * sin_lon) * to_radians / 3600 * 1e9
  
  
  # lod tide
  lod_spline <- interpolate_dut1(self$datetime$utc, 'lod', self$datetime$eop)
  self$lod_t = 1.16 * 2.0 * lod_spline * 
    angular_velocity^2 * earth_radius *
    cos_lat * cos_lat * 1.e9 / 86400.0
  
  
  # earth info
  curvature    <- earth_radius /
    sqrt(1.0 - earth_eccen * sin_lat^2)
  geo_latitude <- 180 / pi * atan(
    ((curvature * (1.0 - earth_eccen) + elevation) * sin_lat) /
      ((curvature + elevation) * cos_lat))                       
  theta <- 90 - geo_latitude
  geo_radius <- sqrt((curvature + elevation)^2 *
                       cos_lat^2 + (curvature * (1.0 - earth_eccen) + elevation)^2 * sin_lat^2)
  df <- 180.0 / pi * 3.600e-3 / gravity
  
  leg <- legendre(6, cos(theta * to_radians))
  
  cos_geo_latitude <- cos(geo_latitude * to_radians)
  sin_geo_latitude <- sin(geo_latitude * to_radians)
  
  radius_ratio <- geo_radius / earth_radius
  
  drdadl <- radius_ratio^(leg[,1])
  dgk <- drdadl * leg[,3]
  dgx <- -1.0 * drdadl / geo_radius * leg[,4] * 1.0e9
  dgy <- drdadl * leg[,2] / (geo_radius * sin(theta * to_radians)) * leg[,3] * 1.0e9
  dgz <- drdadl * leg[,1] / geo_radius * leg[,3] * 1.0e9
  
  dcdlat <- cos_lat * cos_geo_latitude + sin_lat * sin_geo_latitude
  dsdlat <- sin_lat * cos_geo_latitude - cos_lat * sin_geo_latitude
  
  dummy <- dcdlat * dgx - dsdlat * dgz
  dgz <- (dsdlat * dgx + dcdlat * dgz)
  dgx <- dummy
  
  list(azimuth = azimuth,
       latitude = latitude,
       longitude = longitude, 
       elevation = elevation,
       gravity = gravity,
       earth_radius = earth_radius,
       earth_eccen = earth_eccen,
       curvature = curvature,
       geo_latitude = geo_latitude,
       geo_radius = geo_radius,
       theta = theta,
       df = df,
       leg = leg,
       radius_ratio = radius_ratio,
       dgk = dgk, dgx = dgx, 
       dgy = dgy, dgz = dgz,
       angular_velocity = angular_velocity)
}
