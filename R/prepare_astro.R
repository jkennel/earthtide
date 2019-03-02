.prepare_astro <- function(self) {
  
  astro_args <- as.matrix(simon_coef_1994[,-1])
  
  astro <- astro(self$datetime$t_astro,
                 astro_args,
                 self$station$longitude,
                 self$datetime$hours,
                 self$datetime$ddt)
  
  astro_der = astro_der(self$datetime$t_astro, astro_args)
  
  return(list(astro = astro, astro_der = astro_der))
}

