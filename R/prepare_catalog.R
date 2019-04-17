.get_wave_subsets <- function(wave_groups, sub) {
    
  wave_groups <- wave_groups[order(wave_groups$start),]

    freq_sub <- data.frame(
      wave_groups, 
      wave_names = paste0('group_', 
                          sprintf("%02d", 1:nrow(wave_groups))), 
      id = 1:nrow(wave_groups)
    )
    
    # create constituent subgroups
    sub_keep <- list()
    for (i in 1:nrow(freq_sub)) {
      wh <- which(sub$frequency_cpd >= freq_sub$start[i] & 
                    sub$frequency_cpd <= freq_sub$end[i])
      if(length(wh) > 0) {
        
      sub_keep[[i]] <- sub[wh, ]
      sub_keep[[i]]$start <- freq_sub$start[i]
      sub_keep[[i]]$end <- freq_sub$end[i]
      sub_keep[[i]]$wave_names <- freq_sub$wave_names[i]
      sub_keep[[i]]$id <- freq_sub$id[i]
      }
      else {
        sub_keep[[i]] <- NULL
      }
    }
    
    sub <- do.call(rbind, sub_keep)
    
    return(sub)
}

.parse_wave_group <- function(wave_groups = NULL, sub) {
  
  
  # default to range of wave groups
  if (is.null(wave_groups)) {
    return(data.frame(start = 0, end = max(sub$frequency_cpd), multiplier = 1.0))  
  }
  
  
  nc <- ncol(wave_groups)
  # common columns 
  col_matches <- intersect(c('start', 'end', 'multiplier'), names(wave_groups))
  len_matches <- length(col_matches)
  

  if (len_matches == 3L) {
    
    return(wave_groups[, c('start', 'end', 'multiplier')])
    
  } else if (nc == 3L) {
    
    names(wave_groups) <- c('start', 'end', 'multiplier')
    return(wave_groups)
    
  } else if (len_matches == 2L) {
    
    wave_groups$multiplier <- 1.0
    return(wave_groups[, c('start', 'end', 'multiplier')])
    
  } else if (nc == 2L) {
    
    wave_groups$multiplier <- 1.0
    names(wave_groups) <- c('start', 'end', 'multiplier')
    return(wave_groups)
    
  } else if (nc > 3L) {
    
    
    warning('Using the first three columns of data.frame for start, 
             end, and multiplier.  Please either use named columns 
             (start, end, multiplier) or provide a two or three 
             column data.frame if this is not what you desired.')
    
    wave_groups <- wave_groups[, 1:3]
    names(wave_groups) <- c('start', 'end', 'multiplier')
    
    return(wave_groups)
    
  } else {
    
    stop('Wave groups should either be specifed as NULL, or with a data.frame
         with 2 or more columns')
  }
  
}


.prepare_catalog = function(cutoff, wave_groups = NULL, catalog = 'ksm04') {
  
  
  if (catalog == 'hw95s') {
    cat_sub    <- hw95s  
    cat_sub$C2 <- 0.0
    cat_sub$S2 <- 0.0
  } else if (catalog == 'ksm04'){
    cat_sub    <- ksm04
  }
  
  wh  <- which(cat_sub$amplitude > cutoff)
  
  sub <- cat_sub[wh, ]
  
  # combine groups with same astro arguments
  sub <- stats::aggregate(cbind(C0,C1,C2,S0,S1,S2) ~ degree+order+k02+k03+k04+k05+k06+k07+k08+k09+k10+k11+frequency_cpd,
                 data = sub, FUN = sum)
  
  
  wave_groups <- .parse_wave_group(wave_groups, sub)
  sub <- .get_wave_subsets(wave_groups, sub)

  w <- unique(sub[, c('start', 'end')])
  
  
  list(
    catalog = catalog,
    wave_groups = wave_groups,
    k = as.matrix(sub[, c('order', 'k02', 'k03', 'k04',
                              'k05', 'k06', 'k07', 'k08', 'k09',
                              'k10', 'k11')]),
    col_names = paste(c('cos', 'sin'), 
                       rep(paste(w$start, w$end, sep = '_'), each = 2), 
                       sep = '_'),
    cutoff  = cutoff,
    cat_sub = sub,
    degree  = sub$degree,
    order   = sub$order,
    n_constituents = nrow(sub),
    c0 = sub$C0,
    c1 = sub$C1,
    c2 = sub$C2,
    s0 = sub$S0,
    s1 = sub$S1,
    s2 = sub$S2,
    id = as.integer(sub$id),
    jcof = (sub$degree + 1) * sub$degree/2 - 2 + sub$order
  )
  
}