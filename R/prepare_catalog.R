.prepare_catalog = function(cutoff, freq_range, catalog = 'ksm03') {
  
  
  utils::data(list = catalog, package = 'earthtide')
  
  cat_sub <- get(catalog)
  
  if (catalog == 'hw95s') {
    cat_sub$C2 <- 0.0
    cat_sub$S2 <- 0.0
  }
  
  wh  <- which(cat_sub$amplitude > cutoff)
  
  sub <- cat_sub[wh, ]
  
  # combine groups with same astro arguments
  sub <- stats::aggregate(cbind(C0,C1,C2,S0,S1,S2) ~ degree+order+k02+k03+k04+k05+k06+k07+k08+k09+k10+k11+frequency_cpd,
                 data = sub, FUN = sum)
  
  
  if (!all(is.na(freq_range))) {
    freq_range <- freq_range[order(freq_range[, 1]),]
    
    freq_sub <- data.frame(
      freq_range, 
      wave_names = paste0('group_', 
                          sprintf("%02d", 1:nrow(freq_range))), 
      id = 1:nrow(freq_range)
    )
    
    # create constituent subgroups
    sub_keep <- list()
    for (i in 1:nrow(freq_sub)) {
      wh <- which(sub$frequency_cpd >= freq_sub$start[i] & 
                    sub$frequency_cpd <= freq_sub$end[i])
      sub_keep[[i]] <- sub[wh, ]
      sub_keep[[i]]$start <- freq_sub$start[i]
      sub_keep[[i]]$end <- freq_sub$end[i]
      sub_keep[[i]]$wave_names <- freq_sub$wave_names[i]
      sub_keep[[i]]$id <- freq_sub$id[i]
    }
    
    sub <- do.call(rbind, sub_keep)
    
    
  } else {
    sub$id <- 1
    sub$wave_names <- 'group_01'
    sub$start <- 0
    sub$end <- max(sub$frequency_cpd)
    warning('No frequency range specified: defaulting to catalog range')
  }
  
  w <- unique(sub[, c('start', 'end')])
  
  
  list(
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