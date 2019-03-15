# love
#
# returns the love numbers - direct port of eterna 3.40 code
# http://igets.u-strasbg.fr/soft_and_tool.php
#
# @param latitude latitude the station latitude
# @param elevation elevation the station elevation
#
# @return list of love numbers
# 
# @keywords internal
# 
love <- function(latitude, elevation){
  
  dg0 <- c(1.1576,1.1542,1.1600,1.0728,1.0728,1.0728,
           1.0728,1.0363,1.0363,1.0363,1.0363,1.0363)
  dgp <- c(-0.0016,-0.0018,-0.0010,0.0,0.0,0.0,-0.0010,
           0.0,0.0,0.0,0.0,-0.000315)
  dgm <- c(0.0054,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  dh0 <- c(0.6165,0.6069,0.6133,0.2946,0.2946,0.2946,
           0.2946,0.1807,0.1807,0.1807,0.1807,0.1807)
  dhp <- c(0.0007,0.0007,0.0005,0.0,0.0,0.0,0.0003,
           0.0,0.0,0.0,0.0,0.00015)
  dhm <- c(0.0018,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
           0.0,0.0)
  dk0 <- c(0.3068,0.3009,0.3034,0.0942,0.0942,0.0942,
           0.0942,0.0427,0.0427,0.0427,0.0427,0.0427)
  dkp <- c(0.0015,0.0014,0.0009,0.0,0.0,0.0,0.0007,
           0.0,0.0,0.0,0.0,0.00066)
  dkm <- c(-0.0004,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
  
  #Shida-numbers
  dl0   <- c(0.08400,0.08410,0.08520,0.01490,0.01490,0.01490,
             0.01490,0.01000,0.01000,0.01000,0.01000,0.01000)
  dlp   <- c(-0.0020,-0.0020,-0.0010,0.00000,0.00000,0.00000,
             0.00000,0.00000,0.00000,0.00000,0.00000,0.00000)
  dlm   <- c(0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,
             0.00000,0.00000,0.00000,0.00000,0.00000,0.00000)
  latitude_p <- c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
             0.0,0.0)
  latitude_m <- c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
             0.0,0.0)
  
  dea  <- 6378136.3
  dee  <- 6.69439795140e-3
  
  domr <- 15.073729
  dom0 <- 13.943036
  dgr  <- -0.000625
  dhr  <- -0.002505
  dkr  <- -0.001261
  dlr  <- 0.0000781
  
  dclat <- cos(latitude * pi / 180)
  dslat <- sin(latitude * pi / 180)
  dn <- dea / sqrt(1.0 - dee * dslat^2) # curvature
  
  dpsi  <- 180 / pi * atan(((dn*(1.0-dee)+elevation)*dslat)/
                             ((dn+elevation)*dclat)) # geo_latitude
  dthet <- 90.0 - dpsi # theta
  dct   <- cos(dthet * pi / 180)
  dct2  <- dct * dct
  
  latitude_p[1] <- 0.335410 * (35.0 * dct2 * dct2 - 30.0 * dct2 + 3.0) /
                              (3.0 * dct2-1.0)
  
  latitude_m[1]  <- 0.894427 / (3.0 * dct2 - 1.0)
  latitude_p[2]  <- 0.612372 * (7.0 * dct2 - 3.0)
  latitude_p[3]  <- 0.866025 * (7.0 * dct2 - 1.0)
  latitude_p[7]  <- 0.829156 * (9.0 * dct2 - 1.0)
  latitude_p[12] <- 0.806226 * (11.0 * dct2 - 1.0)
  
  dglat <- dg0 + dgp * latitude_p + dgm * latitude_m
  dhlat <- dh0 + dhp * latitude_p + dhm * latitude_m
  dklat <- dk0 + dkp * latitude_p + dkm * latitude_m
  dllat <- dl0 + dlp * latitude_p + dlm * latitude_m
  
  dtlat <- 1.0 + dk0 - dh0 + latitude_p * (dkp-dhp) + latitude_m * (dkm - dhm)
  dtr   <- dkr - dhr      
  
  return(list(dom0=dom0, domr=domr, dgr=dgr, dhr=dhr, dkr=dkr, dlr=dlr, dtr=dtr,
              dglat=dglat, dhlat=dhlat, dklat=dklat, dllat=dllat, dtlat=dtlat))
}
