#' get_main_frequency
#'
#' Get the frequency of the wave with the maximum amplitude in a range.
#'
#' @param start the starting frequency in cycles per day (numeric)
#' @param end the ending frequency in cycles per day (numeric)
#'
#' @return the main frequency between start and end
#' @export
#'
get_main_frequency <- function(start, end) {
  
  if(length(start) != length(end)) {
    stop('start and end must be the same length')  
  }
  
  if(any(start > end)) {
    stop('all values of start must be greater than the corresponding values of
         end')
  }
  
  mn <- c()
  for (i in seq_along(start)) {
    inds <- ksm04$frequency_cpd >= start[i] & ksm04$frequency_cpd < end[i]
    et_sub <- ksm04[inds, ]
    mn[i] <- et_sub[which.max(et_sub$amplitude), 'frequency_cpd']
    
  }
  
  return(mn)
  
}