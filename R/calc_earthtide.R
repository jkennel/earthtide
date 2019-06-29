#' @title earthtide
#' 
#' @description This is a wrapper to the Earthtide R6 class for the prediction 
#' of Earth tides. This function is provided for users who would prefer a more 
#' typical R function.
#'
#' @param utc The date-time in UTC (POSIXct vector).
#' @param do_predict run in predict or analyze mode
#' @param method One or more of "gravity", 
#'     "tidal_potential", "tidal_tilt", "vertical_displacement", 
#'     "horizontal_displacement", "n_s_displacement", "e_w_displacement",
#'     "vertical_strain", "areal_strain", "volume_strain", "horizontal_strain",
#'     or "ocean_tides", "pole_tide", "lod_tide". The pole tide and lod_tide 
#'     are used in predict mode even if do_predict is FALSE. More than one value
#'     can only be used if do_predict == TRUE.
#' @param astro_update Integer that
#'     determines how often to phases are updated in number of samples. Defaults
#'     to 1 (every sample), but speed gains are realized with larger values.
#'     Typically updating every hour will have speed gains and keep precision 
#'     (ie 3600 for one second data, 60 for minute data, 1 for hourly data).
#' @param latitude The station latitude (numeric) defaults to 0.
#' @param longitude The station longitude (numeric) defaults to 0.
#' @param elevation The station elevation (m) (numeric) defaults to 0.
#' @param azimuth Earth azimuth (numeric) defaults to 0.
#' @param gravity Gravity at the station (m/s^2) (numeric) 0 to 
#'     estimate gravity from elevation and latitude.
#' @param earth_radius Radius of earth (m) (numeric) defaults to 6378136.3
#' @param earth_eccen Eccentricity of earth (numeric) 
#'     defaults to 6.69439795140e-3
#' @param cutoff Cutoff amplitude for constituents (numeric) 
#'     defaults to 1e-6.
#' @param wave_groups Two column data.frame having start and end of 
#'     frequency groups (data.frame). This data.frame must have two columns
#'     with the names 'start', and 'end' signifying the start and end of the 
#'     wave groupings.  An optional third column 'multiplier' can be provided 
#'     to scale the particular wave group.  If column names do no match, the
#'     inferred column positions are start, end, multiplier.
#' @param catalog Use the "hw95s" catalog or "ksm04" catalog (character).
#' @param eop User defined Earth Orientation Parameter (EOP) data.frame with the 
#'     following columns: datetime, ddt, ut1_utc, lod, x, y, dx, dy
#' @param return_matrix Return a matrix of tidal values instead of data.frame. 
#'     The datetime column will not be present in this case (logical).
#' @param scale Scale results when do_predict is FALSE
#' @param ... Currently not used.
#'
#' @return data.frame of tidal results
#' @export
#'
#' @examples
#' tms <- as.POSIXct('1990-01-01', tz = 'UTC') + c(0, 3600)
#' wave_groups = data.frame(start = 0, end = 8, multiplier = 1.5)
#' 
#' et <- calc_earthtide(utc = tms,
#'                     do_predict = TRUE,
#'                     method = c('tidal_potential', 'lod_tide', 'pole_tide'),
#'                     astro_update = 1,
#'                     latitude = 52.3868,
#'                     longitude = 9.7144,
#'                     elevation = 110,
#'                     gravity = 9.8127, 
#'                     cutoff = 1.0e-5,
#'                     catalog = 'ksm04',
#'                     wave_groups = wave_groups)
calc_earthtide <- function(utc,
                      do_predict = TRUE,
                      method = 'gravity',
                      astro_update = 1,
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
                      return_matrix = FALSE,
                      scale = TRUE,
                      ...){
  
  
  et <- Earthtide$new(utc = utc, 
                      latitude = latitude,
                      longitude = longitude,
                      elevation = elevation,
                      gravity = gravity, 
                      cutoff = cutoff,
                      catalog = catalog,
                      wave_groups = wave_groups,
                      eop = eop,
                      azimuth = azimuth,
                      earth_radius = earth_radius,
                      earth_eccen = earth_eccen)
  
  if(length(method) > 1) {
    if(!do_predict) {
      stop('If do_predict is FALSE only one method can be provided.')
    }
  }
  
  for (i in seq_along(method)) {
    
    if (method[i] == 'pole_tide') {
      
      et$pole_tide()
      
    } else if (method[i] == 'lod_tide') {
      
      et$lod_tide()
      
    } else if (do_predict) {
      
      if(return_matrix) {
        return(et$predict(method = method[i], 
                          astro_update = astro_update, 
                          return_matrix = return_matrix))  
      } else {
        et$predict(method = method[i], 
                   astro_update = astro_update, 
                   return_matrix = return_matrix)
        
      }
    } else {
      if(return_matrix) {
        return(et$analyze(method = method[i], 
                          astro_update = astro_update, 
                          return_matrix = return_matrix, 
                          scale = scale)
        )
      } else {
        
        et$analyze(method = method[i], 
                   astro_update = astro_update, 
                   return_matrix = return_matrix,
                   scale = scale)
      }
      
    }
    
  }
  
  return(et$tide())
  
}

