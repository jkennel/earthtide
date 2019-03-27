.prepare_datetime <- function(utc, eop = NULL) {
  
  # set the earth orientation parameters
  if (is.data.frame(eop)) {
    if (all(names(dut1) %in% names(eop))){
      eop <- eop
    } else {
      stop(paste0('eop must have the following column names: ', 
                  paste0(names(dut1), collapse = ', ')))
    }
  } else {
    eop <- dut1
  }
  
  utc      <- sort(utc)
  mjd      <- utc_mod_julian(utc)
  j2000    <- utc_julian_2000(utc)
  ddt_int  <- interpolate_dut1(utc, 'ddt', eop) 
  dut1_int <- interpolate_dut1(utc, 'ut1_utc', eop) 
  hours    <- (as.numeric(utc) %% 86400.0 + dut1_int) / 3600.0
  t_astro  <- (j2000 + ddt_int / 3155760000.0) / 10.0

  
  list(utc = utc, mjd = mjd, j2000 = j2000, ddt = ddt_int, dut1 = dut1_int, 
       hours = hours, t_astro = t_astro, eop = eop)
  
}