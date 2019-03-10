#' @title Hartmann and Wenzel (1995) (ETERNA 3.4) wavegroups
#' @format A \code{data.frame} The columns are:
#' \describe{
#'  \item{\code{name}}{wave group name}
#'  \item{\code{start}}{lower frequency of the wave group}
#'  \item{\code{end}}{higher frequency of the wave group}
#'  \item{\code{time}}{applicable to data of what length}
#' }
#' @examples
#' \dontrun{
#' utils::data(eterna_wavegroups)
#' summary(eterna_wavegroups)
#' }
'eterna_wavegroups'

#' @title Hartmann and Wenzel (1995) tidal potential catalogue
#' @format A \code{data.frame} The columns are:
#' \describe{
#'  \item{\code{index}}{the index}
#'  \item{\code{body}}{body generating the potential}
#'  \item{\code{degree}}{degree of the spherical harmonic development}
#'  \item{\code{order}}{order of the spherical harmonic development, k1 = integer argument number for the mean local Moontime, period 24 hours 50 min.}
#'  \item{\code{k02}}{integer argument number for the mean longitude of the Moon, period 27.3 days.}
#'  \item{\code{K03}}{integer argument number for the mean longitude of the Sun, period 365.25 days.}
#'  \item{\code{K04}}{integer argument number for the mean longitude of the lunar perigee, period 8.8 years.}
#'  \item{\code{K05}}{integer argument number for the negate mean longitude of the lunar ascending node, period 18.6 years.}
#'  \item{\code{K06}}{integer argument number for the mean longitude of the solar perigee, period 20942 years.}
#'  \item{\code{K07}}{integer argument number for the mean longitude of Mercury, period 88 days.}
#'  \item{\code{K08}}{integer argument number for the mean longitude of Venus, period 225 days.}
#'  \item{\code{K09}}{integer argument number for the mean longitude of Mars, period 1.88 years.}
#'  \item{\code{K10}}{integer argument number for the mean longitude of Jupiter, period 11.86 years.}
#'  \item{\code{K11}}{integer argument number for the mean longitude of Saturn, period 29.4 years.}
#'  \item{\code{frequency}}{frequency of the tidal wave at J2000 in degree per hour.}
#'  \item{\code{C0}}{COS-coefficient of the tidal potential in 10^-10 m^2/s^2. The C0 coefficient has to be multiplied with the COS of the argument.}
#'  \item{\code{S0}}{SIN-coefficient of the tidal potential in 10^-10 m^2/s^2. The S0 coefficient has to be multiplied with the SIN of the argument.}
#'  \item{\code{C1}}{t*COS-coefficent of the tidal potential in 10^-10 m^2/s^2 per Julian century. The C1 coefficient has to be multiplied with the time difference between the epoch and J2000 (in Julian centuries) and with the COS of the argument.}
#'  \item{\code{S1}}{t*SIN-coefficent of the tidal potential in 10^-10 m^2/s^2 per Julian century. The S1 coefficient has to be multiplied with the time difference between the epoch and J2000 (in Julian centuries) and with the SIN of the argument.}
#'  \item{\code{name}}{Darwin name of the tidal wave (for very few main tidal waves available only)}
#'  \item{\code{frequency_cpd}}{frequency of the tidal wave at J2000 in cycles per day}
#'  \item{\code{amplitude}}{amplitude of the tidal wave at J2000}
#'  \item{\code{phase}}{phase of the tidal wave at J2000}
#' }
#' 
#' @references Hartmann, T., Wenzel, H.-G., 1995. The HW95 tidal potential catalogue. Geophys. Res. Lett. 22, 3553–3556. \doi{10.1029/95GL03324}
#'
#' @examples
#' \donttest{
#' data(hw95s)
#' summary(hw95s)
#' }
'hw95s'




#' @title Kudryavtsev, S M (2004) tidal potential catalogue
#' @format A \code{data.frame} The columns are:
#' \describe{
#'  \item{\code{index}}{the index}
#'  \item{\code{body}}{body generating the potential}
#'  \item{\code{degree}}{degree of the spherical harmonic development}
#'  \item{\code{order}}{order of the spherical harmonic development, k1 = integer argument number for the mean local Moontime, period 24 hours 50 min.}
#'  \item{\code{k02}}{integer argument number for the mean longitude of the Moon, period 27.3 days.}
#'  \item{\code{K03}}{integer argument number for the mean longitude of the Sun, period 365.25 days.}
#'  \item{\code{K04}}{integer argument number for the mean longitude of the lunar perigee, period 8.8 years.}
#'  \item{\code{K05}}{integer argument number for the negate mean longitude of the lunar ascending node, period 18.6 years.}
#'  \item{\code{K06}}{integer argument number for the mean longitude of the solar perigee, period 20942 years.}
#'  \item{\code{K07}}{integer argument number for the mean longitude of Mercury, period 88 days.}
#'  \item{\code{K08}}{integer argument number for the mean longitude of Venus, period 225 days.}
#'  \item{\code{K09}}{integer argument number for the mean longitude of Mars, period 1.88 years.}
#'  \item{\code{K10}}{integer argument number for the mean longitude of Jupiter, period 11.86 years.}
#'  \item{\code{K11}}{integer argument number for the mean longitude of Saturn, period 29.4 years.}
#'  \item{\code{frequency}}{frequency of the tidal wave at J2000 in degree per hour.}
#'  \item{\code{C0}}{COS-coefficient of the tidal potential in 10^-10 m^2/s^2. The C0 coefficient has to be multiplied with the COS of the argument.}
#'  \item{\code{S0}}{SIN-coefficient of the tidal potential in 10^-10 m^2/s^2. The S0 coefficient has to be multiplied with the SIN of the argument.}
#'  \item{\code{C1}}{t*COS-coefficent of the tidal potential in 10^-10 m^2/s^2 per Julian century. The C1 coefficient has to be multiplied with the time difference between the epoch and J2000 (in Julian centuries) and with the COS of the argument.}
#'  \item{\code{S1}}{t*SIN-coefficent of the tidal potential in 10^-10 m^2/s^2 per Julian century. The S1 coefficient has to be multiplied with the time difference between the epoch and J2000 (in Julian centuries) and with the SIN of the argument.}
#'  \item{\code{C2}}{t*COS-coefficent of the tidal potential in 10^-10 m^2/s^2 per Julian century. The C1 coefficient has to be multiplied with the time difference between the epoch and J2000^2 (in Julian centuries) and with the COS of the argument.}
#'  \item{\code{S2}}{t*SIN-coefficent of the tidal potential in 10^-10 m^2/s^2 per Julian century. The S1 coefficient has to be multiplied with the time difference between the epoch and J2000^2 (in Julian centuries) and with the SIN of the argument.}
#'  \item{\code{name}}{Darwin name of the tidal wave (for very few main tidal waves available only)}
#'  \item{\code{amplitude}}{amplitude of the tidal wave at J2000}
#'  \item{\code{frequency_cpd}}{frequency of the tidal wave at J2000 in cycles per day}
#' }
#' 
#' @references Kudryavtsev, S.M., 2004. Improved harmonic development of the Earth tide-generating potential. J. Geod. 77, 829–838. \doi{10.1007/s00190-003-0361-2}
#' 
#' @examples
#' \donttest{
#' utils::data(ksm04)
#' summary(ksm04)
#' }
'ksm04'



#' @title Simon 1994 astronomical constants
#' @format A \code{data.frame} The columns are:
#' \describe{
#'  \item{\code{j}}{7 = mercury, 8 = venus, 9 = mars, 10 = Jupiter, 11 = Saturn}
#'  \item{\code{constant}}{constant}
#'  \item{\code{t1}}{t1}
#'  \item{\code{t2}}{t2}
#'  \item{\code{t3}}{t3}
#'  \item{\code{t4}}{t4}
#' }
#' @references Simon JL, Bretagnon P, Chapront J, Chapront-Touzè M, Francou G, Laskar J (1994) Numerical expressions for precession for- mulae and mean elements for the Moon and planets. Astron Astrophys 282:663–683
#' 
#' @examples
#' \donttest{
#' utils::data(simon_coef_1994)
#' summary(simon_coef_1994)
#' }
'simon_coef_1994'


#' @title dut1 
#' @format A \code{data.frame} The columns are:
#' \describe{
#'  \item{\code{datetime}}{UTC time}
#'  \item{\code{ddt}}{TDT-UTC}
#'  \item{\code{ut1_utc}}{UT1-UTC}
#'  \item{\code{tai_utc}}{TAI - UTC}
#'  \item{\code{lod}}{length of day correction}
#' }
#' 
#' @references http://hpiers.obspm.fr/eop-pc/index.php 
#' 
#' @examples
#' \donttest{
#' utils::data(dut1)
#' }
'dut1'



#' get_dut1
#'
#' Downloads earth orientation data from http://hpiers.obspm.fr/eop-pc/index.php 
#' 
#' @return information on time conversions and leap seconds
#'
get_dut1 <- function(){
 
  tf <- tempfile()
  
  utils::download.file('http://hpiers.obspm.fr/eop-pc/products/combined/C04.php?date=1&eop=7&year1=1962&month1=1&day1=1&year2=2099&month2=12&day2=31&SUBMIT=Submit+Search', tf)
  
  len  <- length(readLines(tf))
  dut1 <- utils::read.table(tf, skip = 2, nrows = len-3, stringsAsFactors = FALSE,
                col.names = c('year', 'month', 'day', 'x', 'x_sig', 'y', 'y_sig',
                              'ut1_utc', 'ut1_utc_sig', 'ut1_tai', 'ut1_tai_sig',
                              'lod', 'lod_sig', 'w3', 'w3_sig',
                              'dx', 'dx_sig', 'dy', 'dy_sig'),
                colClasses = c('integer', 'integer', 'integer',
                               rep('numeric', 16)))
  dut1 <- dut1[-nrow(dut1),]
  
  
  dut1$datetime <- as.POSIXct(paste(dut1$year, 
                                    sprintf("%02i", dut1$month),
                                    sprintf("%02i", dut1$day), sep = '-'),
                              tz = 'UTC')
  
  dut1$ut1_utc <- dut1$ut1_utc / 1000
  dut1$ut1_tai <- dut1$ut1_tai / 1000
  dut1$x       <- dut1$x / 1000
  dut1$y       <- dut1$y / 1000
  dut1$dx      <- dut1$dx / 1000
  dut1$dy      <- dut1$dy / 1000
  dut1$lod     <- dut1$lod / 1000
  
  
  dut1$tai_utc <- dut1$ut1_utc - dut1$ut1_tai
  
  # equation from http://maia.usno.navy.mil/
  # monthly values http://maia.usno.navy.mil/ser7/deltat.data
  dut1$ddt <- 32.184 + (dut1$tai_utc - dut1$ut1_utc)
  
  dut1 <- dut1[, c('datetime', 'ddt', 'ut1_utc', 'tai_utc', 'lod', 'x', 'y',
                   'dx', 'dy')]
  dut1
  
}

# dut1 <- get_dut1()
# usethis::use_data(dut1, internal = FALSE, overwrite = TRUE)
#usethis::use_data(ksm04, internal = FALSE, overwrite = TRUE)
