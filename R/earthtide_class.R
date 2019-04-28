#' @title Earthtide class
#'
#' Class to generate synthetic earthtide signals.
#' 
#' @rdname Earthtide_class
#' 
#' @section Usage:
#' \preformatted{
#' et <- Earthtide$new(
#'   utc = as.POSIXct("2017-01-01", tz = "UTC") + 0:(24 * 7) * 3600, 
#'   latitude = 52.3868,
#'   longitude = 9.7144,
#'   catalog = "ksm04",
#'   wave_groups = data.frame(start = 0.0, end = 6.0))
#' 
#' et$predict(method = "gravity", astro_update = 1)
#' et$analyze(method = "gravity", astro_update = 1)
#' et$lod_tide()
#' et$pole_tide()
#' et$tide()
#' et$print()
#' }
#' 
#' @section Arguments:
#' \code{Earthtide$new}
#' \itemize{
#'   \item{et: }{An \code{Earthtide} object.}
#'   \item{utc: }{The date-time in UTC (POSIXct vector).}
#'   \item{latitude: }{The station latitude (numeric) defaults to 0.}
#'   \item{longitude: }{The station longitude (numeric) defaults to 0.}
#'   \item{elevation: }{The station elevation (m) (numeric) defaults to 0.}
#'   \item{azimuth: }{Earth azimuth (numeric) defaults to 0.}
#'   \item{gravity: }{Gravity at the station (m/s^2) (numeric) 0 to 
#'     estimate gravity from elevation and latitude.}
#'   \item{earth_radius: }{Radius of earth (m) (numeric) defaults to 6378136.3 }
#'   \item{earth_eccen: }{Eccentricity of earth (numeric) 
#'     defaults to 6.69439795140e-3}
#'   \item{cutoff: }{Cutoff amplitude for constituents (numeric) 
#'     defaults to 1e-6}
#'   \item{wave_groups: }{Two column data.frame having start and end of 
#'     frequency groups (data.frame). This data.frame must have two columns
#'     with the names 'start', and 'end' signifying the start and end of the 
#'     wave groupings.  An optional third column 'multiplier' can be provided 
#'     to scale the particular wave group.  If column names do no match, the
#'     inferred column positions are start, end, multiplier.}
#'   \item{catalog: }{Use the "hw95s" catalog or "ksm04" catalog (character).}
#'   \item{eop: }{User defined Earth Orientation Parameter (EOP) data.frame with the 
#'     following columns: datetime, ddt, ut1_utc, lod, x, y, dx, dy}
#'   \item{...: }{Currently not used.}
#' }
#' 
#' \code{Earthtide$predict, Earthtide$analyze}
#' \itemize{
#'   \item{method: }{For \code{predict} and \code{analyze}. One of "gravity", 
#'     "tidal_potential", "tidal_tilt", "vertical_displacement", 
#'     "vertical_strain", "areal_strain", "volume_strain", or "ocean_tides".}
#'   \item{astro_update: }{For \code{predict} and \code{analyze}. Integer that
#'     determines how often to phases are updated in number of samples. Defaults
#'     to 1 (every sample), but speed gains are realized with larger values.
#'     Typically updating every hour will have speed gains and keep precision 
#'     (ie 3600 for one second data, 60 for minute data, 1 for hourly data).}
#' }
#' 
#' @section Details:
#' 
#' \code{$new(utc, latitude, longitude, elevation, azimuth, gravity,} \cr
#' \code{earth_radius, earth_eccen, cutoff, wave_groups, catalog, ...)} \cr
#' create a new \code{Earthtide} object and initialize catalog, station and times.
#' 
#' \code{$predict(method, astro_argument)} generate a combined 
#'   synthetic Earth tide.
#'   
#' \code{$analyze(method, astro_argument)} generate components 
#'   of the Earth tide for analysis.
#'   
#' \code{$lod_tide()} generate components of the LOD (Length Of Day) tide.
#' 
#' \code{$pole_tide()} generate components of the pole tide.
#' 
#' \code{$tide()} get the tide \code{data.frame}.
#' 
#' \code{$print()} print the \code{Earthtide} object.
#' 
#' @docType class
#' @aliases Earthtide-class
#' @importFrom R6 R6Class
#' @importFrom stats approx
#' @importFrom utils read.table
#' @importFrom utils read.fwf
#' @importFrom utils download.file
#' @importFrom utils data
#' @format An \code{\link{R6Class}} generator object
#' 
#' 
#' @references Hartmann, T., Wenzel, H.-G., 1995. The HW95 tidal potential catalogue. Geophys. Res. Lett. 22, 3553-3556. \doi{10.1029/95GL03324}
#' @references Kudryavtsev, S.M., 2004. Improved harmonic development of the Earth tide-generating potential. J. Geod. 77, 829-838. \doi{10.1007/s00190-003-0361-2}
#' @references Wenzel, H.G., 1996. The nanogal software: Earth tide data processing package ETERNA 3.30. Bull. Inf. Marées Terrestres, 124, pp.9425-9439. \url{http://www.eas.slu.edu/GGP/ETERNA34/MANUAL/ETERNA33.HTM}
#' 
#' @examples
#' 
#' et <- Earthtide$new(
#'   utc = as.POSIXct("2017-01-01", tz = "UTC") + 0:(24 * 7) * 3600, 
#'   latitude = 52.3868,
#'   longitude = 9.7144,
#'   catalog = "ksm04",
#'   wave_groups = data.frame(start = 0.0, end = 6.0))
#' 
#' et$predict(method = "gravity", astro_update = 1)
#' 
#' plot(gravity~datetime, et$tide(), type='l')
#' 
#' @name Earthtide
NULL

#' @export

Earthtide <- R6Class("et",
  public = list(

    # initialization
    initialize = function(utc,
                          latitude = 0,
                          longitude = 0, 
                          elevation = 0,
                          azimuth = 0,
                          gravity = 0,
                          earth_radius = 6378136.3,
                          earth_eccen = 6.69439795140e-3,
                          cutoff = 1e-6,
                          wave_groups = NULL,
                          catalog = 'ksm04',
                          eop = NULL,
                          ...) {
      
      # Initialize class using input values
      self$prepare_datetime(utc, eop)
      
      self$prepare_station(latitude, longitude, elevation, azimuth, gravity,
                           earth_radius, earth_eccen)

      self$prepare_astro()

      self$prepare_catalog(cutoff, wave_groups, catalog)

      self$love_params <- love(latitude, elevation)
      
      invisible(self)
    },
    
    # Initialize class using input values
    prepare_datetime = function(utc, eop) {
      self$datetime <- .prepare_datetime(utc, eop)
      self$tides <- data.frame(datetime = utc)
    },
    
    # Calculate the astronical arguments and derivative
    prepare_astro = function() {
      self$astro <- .prepare_astro(self)
    },
    
    # Subset values based using a cutoff amplitude and 
    # cutoff frequencies
    prepare_catalog = function(cutoff, wave_groups = NULL, catalog = 'ksm04') {
      
      self$catalog <- .prepare_catalog(cutoff, wave_groups, catalog = catalog)
      
    },
    
    # Calculate properties based on the location of station
    prepare_station = function(latitude, longitude, elevation, azimuth, 
                               gravity, earth_radius, earth_eccen) {
      
      self$station <- .prepare_station(self, latitude, longitude, elevation,
                                       azimuth, gravity, earth_radius,
                                       earth_eccen)
      self$pole_quotient = 1.16 #1.1788, # Ducarme et al. 2006
      self$pk     = rep(0, 25)
      self$delta  = rep(1.0, 25)
      self$deltar = 0.0
      
    }, 
    check_time_increment = function(astro_update) { 
      
      dt <- unique(diff(as.numeric(self$datetime$utc)))
      
      if(length(self$datetime$utc) == 1) {
        self$update_coef <- 0.0
        return(1L)
      }
      
      if (length(dt) != 1L) {
        if (astro_update != 1L) {
          warning('Times are not regularly spaced, setting astro_update to 1') 
        }
        astro_update <- 1L
        self$update_coef <- 0.0
      } else {
        self$update_coef <- pi / 180.0 * dt / 3600.0
        astro_update <- astro_update
      }
      
      astro_update
      
    },
    gravity = function(){

      self$station$dgk <- self$station$dgz
      self$pk[] <- 180 
      
      self$delta[1:12] <- self$love_params$dglat
      self$deltar <- self$love_params$dgr
      
    },
    tidal_potential = function() {
      
      self$delta[1:12] <- self$love_params$dklat
      self$deltar <- self$love_params$dkr

    },
    tidal_tilt = function() {
      
      cos_azimuth <- cos(self$station$azimuth)
      sin_azimuth <- sin(self$station$azimuth)
      x_comp <- self$station$dgx[1:12] * cos_azimuth
      y_comp <- self$station$dgy[1:12] * sin_azimuth
      self$station$dgk[1:12] <- sqrt((x_comp)^2 + (y_comp)^2) * self$station$df
      wh <- which(x_comp != 0 | y_comp != 0)
      self$pk[wh] <- 180 / pi * atan2(y_comp, x_comp)
      
      # from etpots
      self$delta[1:12] <- self$love_params$dtlat
      self$deltar = self$love_params$dkr - self$love_params$dhr
      
    }, 
    vertical_displacement = function() {
      dfak <- 1e3 / self$station$gravity
      self$station$dgk[1:12] <- self$station$dgk[1:12] * 
        self$love_params$dhlat[1:12] * dfak
      self$pk[] <- 0.0
    },
    # This number is way too big from eterna - currently must be an error
    # horizontal_displacement = function() {
    #   
    #   #dfak <- 1e3 *  self$station$geo_radius / self$station$gravity
    #   
    #   dfak <- 1
    #   cos_azimuth <- cos(self$station$azimuth)
    #   sin_azimuth <- sin(self$station$azimuth)
    #   x_comp <- self$station$dgx[1:12] * cos_azimuth
    #   y_comp <- self$station$dgy[1:12] * sin_azimuth
    #   self$station$dgk[1:12] <- sqrt((x_comp)^2 + (y_comp)^2) * 
    #     self$love_params$dllat[1:12] * dfak
    #   self$pk[] <- 0.0
    #   wh <- which(x_comp != 0 | y_comp != 0)
    #   self$pk[wh] <- 180 / pi * atan2(y_comp, x_comp)
    #   
    # },
    vertical_strain = function(poisson = 0.25) {
      dfak <- 1.e9 * poisson / (poisson - 1.0)
      self$strain(dfak)
    },
    areal_strain = function() {
      dfak <- 1.e9
      self$strain(dfak)
    },
    volume_strain = function(poisson = 0.25) {
      dfak <- 1.e9 * (1.0 - 2.0 * poisson) / (1.0 - poisson)
      self$strain(dfak)
    },
    strain = function(dfak) {
      
      scale <- c(rep(-6, 3), rep(-12, 4), rep(-20, 5))
      
      self$station$dgk[1:12] <- self$station$dgk[1:12] * dfak *
        (2.0 * self$love_params$dhlat[1:12] + scale * self$love_params$dllat[1:12]) / 
        (self$station$gravity * self$station$geo_radius)
      
    },
    ocean_tides = function() {
      dfak <- 1e3 / self$station$gravity
      self$station$dgk <- self$station$dgk * dfak
      self$pk[] <- 0.0
    },
    predict = function(method = 'gravity', astro_update = 1L, return_matrix = FALSE) {
      
      self$apply_method(method) 
      astro_update <- self$check_time_increment(astro_update)
      if (return_matrix) {
        mat <- self$calculate(astro_update = astro_update, predict = TRUE)
        colnames(mat) <- method
        return(mat)
      } else {
        self$tides[[method]] <- 
          as.numeric(self$calculate(astro_update = astro_update, predict = TRUE))
      }
      
      # reset parameters after calculation
      self$prepare_station(self$station$latitude, 
                           self$station$longitude, 
                           self$station$elevation,
                           self$station$azimuth, 
                           self$station$gravity, 
                           self$station$earth_radius,
                           self$station$earth_eccen)
      invisible(self)
    },
    analyze = function(method = 'gravity', astro_update = 1L, return_matrix = FALSE) {
      
      self$apply_method(method)
      astro_update <- self$check_time_increment(astro_update)
      
      if (return_matrix) {
        mat <- self$calculate(astro_update = astro_update, predict = FALSE)
        colnames(mat) <- self$catalog$col_names
        return(mat)
      } else {
        self$tides[self$catalog$col_names] <- 
          self$calculate(astro_update = astro_update, predict = FALSE)
      }
      
      # reset parameters after calculation
      self$prepare_station(self$station$latitude, 
                           self$station$longitude, 
                           self$station$elevation,
                           self$station$azimuth, 
                           self$station$gravity, 
                           self$station$earth_radius,
                           self$station$earth_eccen)
      invisible(self)
    },
    apply_method = function(method) {
      
      if (method == 'gravity') {
        self$gravity()
      } else if (method == 'ocean_tides') {
        self$ocean_tides()
      } else if (method == 'tidal_potential') {
        self$tidal_potential()
      } else if (method == 'tidal_tilt') {
        self$tidal_tilt()
      } else if (method == 'vertical_displacement') {
        self$vertical_displacement()
      # } else if (method == 'horizontal_displacement') {
      #   self$horizontal_displacement()
      } else if (method == 'vertical_strain') {
        self$vertical_strain()
      } else if (method == 'areal_strain') {
        self$areal_strain()
      } else if (method == 'volume_strain') {
        self$volume_strain()
      } 
    },
    calculate = function(astro_update = 1L, predict = TRUE) {
      et_calculate(self$astro$astro,
                 self$astro$astro_der,
                 self$catalog$k,
                 self$pk,
                 self$delta,
                 self$deltar,
                 self$catalog$c0,
                 self$catalog$s0,
                 self$catalog$c1,
                 self$catalog$s1,
                 self$catalog$c2,
                 self$catalog$s2,
                 self$station$dgk,
                 self$catalog$jcof-1,
                 self$datetime$j2000,
                 self$love_params$dom0,
                 self$love_params$domr,
                 self$catalog$id,
                 astro_update,
                 self$update_coef, 
                 self$catalog$wave_groups$multiplier,
                 predict)
    },
    pole_tide = function() {
      self$tides$pole_tide <- self$pole_t
      invisible(self)
    },
    lod_tide = function() {
      self$tides$lod_tide <- self$lod_t
      invisible(self)
    },
    tide = function() {
      invisible(self$tides)
    },
    print = function(...) {
      cat("Earthtide: \n")
      cat("  Times: ", length(self$datetime$utc), "\n", sep = '')
      cat("    range (utc): ", as.character(head(self$datetime$utc, 1)), ' to ', 
                               as.character(tail(self$datetime$utc, 1)), "\n", sep = "")
      cat("    dt[1] (sec): ", as.numeric(self$datetime$utc[2]) -
                               as.numeric(self$datetime$utc[1]), "\n", sep = '')
      cat("  Station: \n")
      cat("    latitude:    ", head(self$station$latitude), "\n", sep = "")
      cat("    longitude:   ", head(self$station$longitude), "\n", sep = "")
      cat("    elevation:   ", head(self$station$elevation), "\n", sep = "")
      cat("    gravity:     ", head(self$station$gravity), "\n", sep = "")
      cat("    azimuth:     ", head(self$station$azimuth), "\n", sep = "")
      cat("  Wave groups: \n")
      cat("    catalog:     ", head(self$catalog$catalog), "\n", sep = "")
      cat("    cutoff:      ", head(self$catalog$cutoff), "\n", sep = "")
      cat("    n waves:     ", head(self$catalog$n_constituents), "\n", sep = "")
      cat("    n groups:    ", nrow(self$catalog$wave_groups), "\n", sep = "")
      # cat("  Tides: \n")
      # print(head(self$tides), row.names = FALSE)
      
      invisible(self)
    },
    
    # time variables
    datetime = list(),
    update_coef = NA_real_, 
    
    # astro arguments
    astro = list(),
    
    # wavegroup catalog
    catalog = list(),
    
    # geodetic coeficients
    station = list(),
    
    # love and shida numbers
    love_params = list(),
    
    pole_quotient = 1.16, #1.1788, # Ducarme et al. 2006
    pk     = rep(0, 25),
    delta  = rep(1.0, 25),
    deltar = 0.0,
    
    # outputs
    lod_t = NA_real_,
    pole_t = NA_real_,
    tides = list()
    
  )
)




