.prepare_datetime <- function(utc) {
  
  utc     <- sort(utc)
  mjd     <- utc_mod_julian(utc)
  j2000   <- utc_julian_2000(utc)
  ddt     <- interpolate_dut1(utc, 'ddt') 
  dut1    <- interpolate_dut1(utc, 'ut1_utc') 
  hours   <- (as.numeric(utc) %% 86400.0 + dut1) / 3600.0
  t_astro <- (j2000 + ddt / 3155760000.0) / 10.0

  
  list(utc = utc, mjd = mjd, j2000 = j2000, ddt = ddt, dut1 = dut1, 
       hours = hours, t_astro = t_astro)
  
}