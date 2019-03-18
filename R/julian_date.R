# utc_julian
#
# \code{utc_julian} returns the julian datetime from UTC
#
# @param utc datetime in utc
#
# @return julian datetime
# 
# @keywords internal
# 
utc_julian <- function(utc) {
  
  (as.numeric(utc) / 86400.0) + 2440587.5;
  
}


# utc_julian_2000
#
# \code{utc_julian_2000} returns the j2000 datetime from UTC
#
# @param utc datetime in utc
#
# @return julian datetime
# 
# @keywords internal
# 
utc_julian_2000 <- function(utc) {
  
  (utc_julian(utc) - 2451545.0) / 36525.0

}

# utc_mod_julian
#
# \code{utc_mod_julian} returns the modified julian datetime from UTC
#
# @param utc datetime in utc
#
# @return julian datetime
# 
# @keywords internal
# 
utc_mod_julian <- function(utc) {
  
  (utc_julian(utc) - 2400000.5)
  
}


# julian_mod_julian
#
# \code{julian_mod_julian} returns the modified julian date from the julian date
#
# @param julian datetime in julian
#
# @return modified julian date
# 
# @keywords internal
# 
julian_mod_julian <- function(julian) {
  
  (julian - 2400000.5)
  
}


# mod_julian_julian
#
# \code{mod_julian_julian} returns the julian date from the modified julian date
#
# @param mod_julian datetime in mod_julian
#
# @return julian date
# 
# @keywords internal
# 
mod_julian_julian <- function(mod_julian) {
  
  (mod_julian + 2400000.5)
  
}

# mod_julian_utc
#
# \code{mod_julian_utc} returns the UTC from the modified julian 
#
# @param mod_julian datetime in mod_julian
#
# @return UTC datetime
# 
# @keywords internal
# 
mod_julian_utc <- function(mod_julian) {

  # (-2440587.5 + 2400000.5) = -40587
  as.POSIXct((mod_julian - 40587) * 86400.0, tz = 'UTC', origin = '1970-01-01')

}

