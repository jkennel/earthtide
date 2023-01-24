#' @title Hartmann and Wenzel (1995) (ETERNA 3.4) wavegroups
#'
#' @description This data.frame contains wavegroups for different data time
#' spans.  The wavegroups should be subset prior to use and the 'time' column
#' provides guidelines based on your input time span.
#'
#' @format A \code{data.frame} The columns are:
#' \describe{
#'  \item{\code{name}}{wave group name}
#'  \item{\code{start}}{lowest frequency of the wave group}
#'  \item{\code{end}}{highest frequency of the wave group}
#'  \item{\code{time}}{applicable to data of what length}
#' }
#'
#' @examples
#' utils::data(eterna_wavegroups)
"eterna_wavegroups"

# @title Hartmann and Wenzel (1995) tidal potential catalogue
# @format A \code{data.frame} The columns are:
# \describe{
#  \item{\code{index}}{the index}
#  \item{\code{body}}{body generating the potential}
#  \item{\code{degree}}{degree of the spherical harmonic development}
#  \item{\code{order}}{order of the spherical harmonic development, k1 = integer argument number for the mean local Moontime, period 24 hours 50 min.}
#  \item{\code{k02}}{integer argument number for the mean longitude of the Moon, period 27.3 days.}
#  \item{\code{K03}}{integer argument number for the mean longitude of the Sun, period 365.25 days.}
#  \item{\code{K04}}{integer argument number for the mean longitude of the lunar perigee, period 8.8 years.}
#  \item{\code{K05}}{integer argument number for the negate mean longitude of the lunar ascending node, period 18.6 years.}
#  \item{\code{K06}}{integer argument number for the mean longitude of the solar perigee, period 20942 years.}
#  \item{\code{K07}}{integer argument number for the mean longitude of Mercury, period 88 days.}
#  \item{\code{K08}}{integer argument number for the mean longitude of Venus, period 225 days.}
#  \item{\code{K09}}{integer argument number for the mean longitude of Mars, period 1.88 years.}
#  \item{\code{K10}}{integer argument number for the mean longitude of Jupiter, period 11.86 years.}
#  \item{\code{K11}}{integer argument number for the mean longitude of Saturn, period 29.4 years.}
#  \item{\code{frequency}}{frequency of the tidal wave at J2000 in degree per hour.}
#  \item{\code{C0}}{COS-coefficient of the tidal potential in 10^-10 m^2/s^2. The C0 coefficient has to be multiplied with the COS of the argument.}
#  \item{\code{S0}}{SIN-coefficient of the tidal potential in 10^-10 m^2/s^2. The S0 coefficient has to be multiplied with the SIN of the argument.}
#  \item{\code{C1}}{t*COS-coefficent of the tidal potential in 10^-10 m^2/s^2 per Julian century. The C1 coefficient has to be multiplied with the time difference between the epoch and J2000 (in Julian centuries) and with the COS of the argument.}
#  \item{\code{S1}}{t*SIN-coefficent of the tidal potential in 10^-10 m^2/s^2 per Julian century. The S1 coefficient has to be multiplied with the time difference between the epoch and J2000 (in Julian centuries) and with the SIN of the argument.}
#  \item{\code{name}}{Darwin name of the tidal wave (for very few main tidal waves available only)}
#  \item{\code{frequency_cpd}}{frequency of the tidal wave at J2000 in cycles per day}
#  \item{\code{amplitude}}{amplitude of the tidal wave at J2000}
#  \item{\code{phase}}{phase of the tidal wave at J2000}
# }
#
# @keywords internal
#
# @references Hartmann, T., Wenzel, H.-G., 1995. The HW95 tidal potential catalogue. Geophys. Res. Lett. 22, 3553-3556. \doi{10.1029/95GL03324}
#
# @examples
# hw95s
# 'hw95s'



# @title Kudryavtsev, S M (2004) tidal potential catalogue
# @format A \code{data.frame} The columns are:
# \describe{
#  \item{\code{index}}{the index}
#  \item{\code{body}}{body generating the potential}
#  \item{\code{degree}}{degree of the spherical harmonic development}
#  \item{\code{order}}{order of the spherical harmonic development, k1 = integer argument number for the mean local Moontime, period 24 hours 50 min.}
#  \item{\code{k02}}{integer argument number for the mean longitude of the Moon, period 27.3 days.}
#  \item{\code{K03}}{integer argument number for the mean longitude of the Sun, period 365.25 days.}
#  \item{\code{K04}}{integer argument number for the mean longitude of the lunar perigee, period 8.8 years.}
#  \item{\code{K05}}{integer argument number for the negate mean longitude of the lunar ascending node, period 18.6 years.}
#  \item{\code{K06}}{integer argument number for the mean longitude of the solar perigee, period 20942 years.}
#  \item{\code{K07}}{integer argument number for the mean longitude of Mercury, period 88 days.}
#  \item{\code{K08}}{integer argument number for the mean longitude of Venus, period 225 days.}
#  \item{\code{K09}}{integer argument number for the mean longitude of Mars, period 1.88 years.}
#  \item{\code{K10}}{integer argument number for the mean longitude of Jupiter, period 11.86 years.}
#  \item{\code{K11}}{integer argument number for the mean longitude of Saturn, period 29.4 years.}
#  \item{\code{frequency}}{frequency of the tidal wave at J2000 in degree per hour.}
#  \item{\code{C0}}{COS-coefficient of the tidal potential in 10^-10 m^2/s^2. The C0 coefficient has to be multiplied with the COS of the argument.}
#  \item{\code{S0}}{SIN-coefficient of the tidal potential in 10^-10 m^2/s^2. The S0 coefficient has to be multiplied with the SIN of the argument.}
#  \item{\code{C1}}{t*COS-coefficent of the tidal potential in 10^-10 m^2/s^2 per Julian century. The C1 coefficient has to be multiplied with the time difference between the epoch and J2000 (in Julian centuries) and with the COS of the argument.}
#  \item{\code{S1}}{t*SIN-coefficent of the tidal potential in 10^-10 m^2/s^2 per Julian century. The S1 coefficient has to be multiplied with the time difference between the epoch and J2000 (in Julian centuries) and with the SIN of the argument.}
#  \item{\code{C2}}{t*COS-coefficent of the tidal potential in 10^-10 m^2/s^2 per Julian century. The C1 coefficient has to be multiplied with the time difference between the epoch and J2000^2 (in Julian centuries) and with the COS of the argument.}
#  \item{\code{S2}}{t*SIN-coefficent of the tidal potential in 10^-10 m^2/s^2 per Julian century. The S1 coefficient has to be multiplied with the time difference between the epoch and J2000^2 (in Julian centuries) and with the SIN of the argument.}
#  \item{\code{name}}{Darwin name of the tidal wave (for very few main tidal waves available only)}
#  \item{\code{amplitude}}{amplitude of the tidal wave at J2000}
#  \item{\code{frequency_cpd}}{frequency of the tidal wave at J2000 in cycles per day}
# }
#
# @references Kudryavtsev, S.M., 2004. Improved harmonic development of the Earth tide-generating potential. J. Geod. 77, 829-838. \doi{10.1007/s00190-003-0361-2}
#
# @keywords internal
#
# @examples
# ksm04
# 'ksm04'



# @title Simon 1994 astronomical constants
# @format A \code{data.frame} The columns are:
# \describe{
#  \item{\code{j}}{7 = mercury, 8 = venus, 9 = mars, 10 = Jupiter, 11 = Saturn}
#  \item{\code{constant}}{constant}
#  \item{\code{t1}}{t1}
#  \item{\code{t2}}{t2}
#  \item{\code{t3}}{t3}
#  \item{\code{t4}}{t4}
# }
# @references Simon JL, Bretagnon P, Chapront J, Chapront-Touzè M, Francou G, Laskar J (1994) Numerical expressions for precession formulae and mean elements for the Moon and planets. Astron Astrophys 282:663-683
#
# @keywords internal
#
# @examples
# simon_coef_1994
#
# 'simon_coef_1994'


# @title dut1
# @format A \code{data.frame} The columns are:
# \describe{
#  \item{\code{datetime}}{UTC time}
#  \item{\code{ddt}}{TDT-UTC}
#  \item{\code{ut1_utc}}{UT1-UTC}
#  \item{\code{tai_utc}}{TAI - UTC}
#  \item{\code{lod}}{length of day correction}
# }
#
# @references http://hpiers.obspm.fr/eop-pc/index.php
#
# @keywords internal
#
# @examples
# dut1
# 'dut1'



# download leap second data
get_tai_utc <- function(tai_utc_path) {
  tf <- tempfile()
  utils::download.file(tai_utc_path, tf) # this may need to be updated

  widths <- c(17, 9, 10, 12, 12, 6, 4, 9, 1)
  tai_utc <- read.fwf(tf, widths = widths, stringsAsFactors = FALSE)
  tai_utc <- tai_utc[, c(2, 4, 6, 8)]
  names(tai_utc) <- c("jd", "tai_utc", "minus_date", "factor")

  tai_utc$mjd <- julian_mod_julian(tai_utc$jd)

  tai_utc
}

# get tai - utc
mjd_tai_utc <- function(mjd, tai_utc_path) {
  tai_utc <- get_tai_utc(tai_utc_path)

  tu <- c()

  for (i in seq_along(mjd)) {
    di <- mjd[i] - tai_utc$mjd

    wh <- max(which(di >= 0))

    tu[i] <- tai_utc$tai_utc[wh] + (mjd[i] - tai_utc$minus_date[wh]) *
      tai_utc$factor[wh]
  }

  tu
}


# Bulletin B
get_iers_b <- function(b_path, tai_utc_path) {
  tf <- tempfile()
  utils::download.file(b_path, tf)

  dut1 <- utils::read.table(tf,
    skip = 14, stringsAsFactors = FALSE,
    col.names = c(
      "year", "month", "day", "mjd",
      "x", "y", "ut1_utc", "lod",
      "dx", "dy", "x_sig", "y_sig",
      "ut1_utc_sig", "lod_sig",
      "dx_sig", "dy_sig"
    ),
    colClasses = c(
      "integer", "integer", "integer",
      rep("numeric", 13)
    )
  )

  dut1$datetime <- mod_julian_utc(dut1$mjd)

  # equation from http://maia.usno.navy.mil/
  dut1$ddt <- 32.184 + (mjd_tai_utc(dut1$mjd, tai_utc_path) - dut1$ut1_utc)


  dut1 <- dut1[, c(
    "datetime", "ddt", "ut1_utc", "lod", "x", "y",
    "dx", "dy"
  )]
}


# Bulletin A
get_iers_a <- function(a_path, daily_path, tai_utc_path) {
  widths <- c(
    2, 2, 2, 1, 8, 1, 1, 1, 9, 9, 1, 9, 9, 2, 1,
    10, 10, 1, 7, 7, 2, 1, 1, 9, 9, 1, 9, 9, 10,
    10, 11, 10, 10
  )

  # historical
  tf_all <- tempfile()
  utils::download.file(a_path, tf_all)
  iers_all <- read.fwf(tf_all, widths = widths, stringsAsFactors = FALSE)
  wh <- c(5, 9, 12, 16, 19, 24, 27)
  wh_names <- c("mjd", "x", "y", "ut1_utc", "lod", "dx", "dy")
  iers_all <- iers_all[, wh]
  names(iers_all) <- wh_names
  iers_all$datetime <- mod_julian_utc(iers_all$mjd)
  iers_all$dx <- iers_all$dx / 1000.0
  iers_all$dy <- iers_all$dy / 1000.0
  iers_all$lod <- iers_all$lod / 1000.0
  iers_all$ddt <- 32.184 + (mjd_tai_utc(iers_all$mjd, tai_utc_path) -
    iers_all$ut1_utc)
  iers_all <- iers_all[, c(
    "datetime", "ddt", "ut1_utc", "lod", "x", "y",
    "dx", "dy"
  )]


  # daily set for update
  tf_daily <- tempfile()
  utils::download.file(daily_path, tf_daily)
  iers_daily <- read.fwf(tf_daily, widths = widths, stringsAsFactors = FALSE)
  iers_daily <- iers_daily[, wh]
  names(iers_daily) <- wh_names
  iers_daily$datetime <- mod_julian_utc(iers_daily$mjd)
  iers_daily$dx <- iers_daily$dx / 1000.0
  iers_daily$dy <- iers_daily$dy / 1000.0
  iers_daily$lod <- iers_daily$lod / 1000.0
  iers_daily$ddt <- 32.184 + (mjd_tai_utc(iers_daily$mjd, tai_utc_path) -
    iers_daily$ut1_utc)
  iers_daily <- iers_daily[, c(
    "datetime", "ddt", "ut1_utc", "lod", "x", "y",
    "dx", "dy"
  )]

  iers_all <- iers_all[iers_all$datetime < min(iers_daily$datetime), ]
  return(rbind(iers_all, iers_daily))
}


#

#' get_iers
#'
#' \code{get_iers} returns a \code{data.frame} of earth orientation
#' parameters from (1962-present).  This function requires an active internet
#' connection. Bulletins A and B are combined giving precedence to B.
#' Approximately (~ 7 MB) of data are downloaded. This function is brittle and
#' may fail when data sources change.
#'
#' @param a_path ftp or http path to download IERS bullitin A
#' @param b_path ftp or http path to download IERS bullitin B
#' @param daily_path ftp or http path to download IERS daily data
#' @param tai_utc_path ftp or http path to tai-utc data
#'
#' @return \code{data.frame} of earth orientation parameters with the following
#' columns: datetime, ddt, ut1_utc, lod, x, y, dx, dy.
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' eop <- get_iers()
#' }
#'
get_iers <- function(
    a_path = NULL,
    b_path = NULL,
    daily_path = NULL,
    tai_utc_path = NULL) {
  if (is.null(a_path)) {
    a_path <- "https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all"
  }
  if (is.null(b_path)) {
    b_path <- "http://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now"
  }
  if (is.null(daily_path)) {
    daily_path <- "https://datacenter.iers.org/data/latestVersion/10_FINALS.DATA_IAU2000_V2013_0110.txt"
  }
  if (is.null(tai_utc_path)) {
    tai_utc_path <- "http://astroutils.astronomy.ohio-state.edu/time/tai-utc.txt"
  }
  # http://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now
  # ftp://ftp.iers.org/products/eop/rapid/standard/finals2000A.all
  # ftp://ftp.iers.org/products/eop/rapid/daily/finals2000A.daily

  bull_a <- get_iers_a(a_path, daily_path, tai_utc_path) # bulletin A
  bull_b <- get_iers_b(b_path, tai_utc_path) # bulletin B

  bull_a <- bull_a[bull_a$datetime > max(bull_b$datetime), ]
  bull_ab <- rbind(bull_b, bull_a)
  bull_ab <- bull_ab[rowSums(is.na(bull_ab[, 2:8])) < 7, ]
  bull_ab[is.na(bull_ab)] <- 0
  bull_ab
}




# library(earthtide)
# # a_path = 'https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all'
# # b_path ='http://hpiers.obspm.fr/iers/eop/eopc04/eopc04_IAU2000.62-now'
# # daily_path = 'https://datacenter.iers.org/data/latestVersion/10_FINALS.DATA_IAU2000_V2013_0110.txt'
# # tai_utc_path = 'http://astroutils.astronomy.ohio-state.edu/time/tai-utc.txt'
# dut1 <- get_iers()
# simon_coef_1994 <- earthtide:::simon_coef_1994
# ksm04 <- earthtide:::ksm04
# hw95s <- earthtide:::hw95s
# usethis::use_data(dut1,
#                   simon_coef_1994,
#                   ksm04,
#                   hw95s,
#                   internal = TRUE,
#                   overwrite = TRUE)
