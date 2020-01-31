library(Rcpp)
library(RcppArmadillo)
sourceCpp("bee_movement.cpp")

## Main function 
#
# n_bees        : numeric vector depicting the different abundances to be assessed
# plot_size     : the number of rows and columns (one number for a square matrix) of the mosaic
# steps         : the number of steps of an individual simulation, equivalent to the sampling time at the field
# mosaic_reps   : number of repetitions of a movement simulation on the same mosaic
# mode          : movement pattern, "random" or "nonrandom"
# mosaic_values : if mode = "nonrandom", a numeric vector with the possible mosaic values that may appear
# probs         : if mode = "nonrandom", a vector of probabilities of occurence of the mosaic values
goforit <- function(n_bees = 5, plot_size = 10, steps = 100, mosaic_reps = 5, mode = "random", mosaic_values = 1:10, 
                     probs = sort(seq(0.025, 0.25, 0.025))){
  
  # Define objects
  total_sampled <- NULL
  first_t <- NULL

  # Define objects for 'hotspots': positions in matrix where bees move most, 
  # for three reference abundances (mean, median and max)
  m_hspots <- matrix(0, nrow = plot_size, ncol = plot_size)
  hotspots <- list(min_abundance    = m_hspots,
                   median_abundance = m_hspots,
                   max_abundance    = m_hspots)
  # Define object for saving the mosaics
  list_mosaics <- vector("list", length(n_bees))
  order <- 0
  # Define object for saving the path of individual bee
  paths <- vector("list", length(mosaic_reps))
  foraging_area <- rep(NA, mosaic_reps)

  pb <- txtProgressBar(max = length(n_bees)*mosaic_reps*steps, style = 3, char = "o", width = 50)
  # Run simulation for each bee abundance
  for(b in 1:length(n_bees)){
    
    if(b %in% c(min(1:length(n_bees)), round(median(1:length(n_bees))), max(1:length(n_bees)))){
      order <- order + 1
    }

    target <- rep(ceiling(plot_size/2), 2)
    # Generate landscape of movement probabilities
    # random = equal probabilities
    # non random = unequeal probabilities
    if(mode == "random"){
      mosaic <- matrix(1, nrow = plot_size, ncol = plot_size)
    } else {
      mosaic <- matrix(sample(mosaic_values, size = plot_size*plot_size, replace = T, prob = probs),
                       nrow = plot_size, ncol = plot_size)
      
      mosaic[target[1] + 1, target[2] + 1] <- max(mosaic_values)
    }
    # Assign a zero probability to the cells located at the edges of the mosaic
    mosaic <- cbind(0, rbind(0, mosaic, 0), 0)
    list_mosaics[[b]]        <- mosaic[2:plot_size, 2:plot_size]
    names(list_mosaics)[[b]] <- n_bees[b]

    # Define objects
    sampled   <- rep(NA, mosaic_reps)
    approachs <- rep(NA, mosaic_reps)
    # Run simulation with the same mosaic (equivalent to number of repetitions of a particular plot at the field)
    for(r in 1:mosaic_reps){
      b_pos <- sample(1:(plot_size*plot_size), n_bees[b], replace =  F)
      bees <- matrix(c(ceiling(b_pos/plot_size), b_pos-plot_size*(ceiling(b_pos/plot_size)-1)), ncol = 2)

      n_sampled <- 0
      first_approach <- NaN
      ibees <- -1
      
      # Create a data.frame containing the position of an individual bee
      if(b == 1) paths[[r]] <- data.frame(x = rep(NA, steps), y = rep(NA, steps))

      # Individual simulation, with a number of steps (equivalent to sampling time of an individual plot at the field)
      for(t in 1:steps){
        # bee_movement: C++ code
        bee_sim <- bee_movement(bees, mosaic, ibees, target, n_sampled, first_approach, t) 
        bees           <- bee_sim[[1]]
        n_sampled      <- bee_sim[[2]]
        first_approach <- bee_sim[[3]]
        ibees          <- bee_sim[[4]]
        # Hotspots
        if(order %in% 1:3){
          for(i in 1:nrow(bees)){
            hotspots[[order]][bees[i, 1], bees[i, 2]] <- hotspots[[order]][bees[i, 1], bees[i, 2]] + 1
          }
        }
        # Path recording
        if(b == 1){
          paths[[r]][t, 1] <- bees[1, 1]
          paths[[r]][t, 2] <- bees[1, 2]
        }
        
        setTxtProgressBar(pb, pb$getVal() + 1)
      }
      
      if(b == 1){
        path_matrix <- matrix(0, ncol = plot_size, nrow = plot_size)
        for(i in 1:nrow(paths[[r]])){
          path_matrix[paths[[r]][i, 1], paths[[r]][i, 2]] <- path_matrix[paths[[r]][i, 1], paths[[r]][i, 2]] + 1
        }
        foraging_area[r] <- (max(paths[[r]]$y) - min(paths[[r]]$y)) * (max(paths[[r]]$x) - min(paths[[r]]$x))
        paths[[r]] <- path_matrix
      }
      
      sampled[r]   <- n_sampled
      approachs[r] <- first_approach
    }
    total_sampled <- c(total_sampled, sampled)
    first_t <- c(first_t, approachs)
  }
  df <- data.frame(total_sampled = total_sampled, 
                   fid           = first_t,
                   abundance     = rep(n_bees, each = mosaic_reps))
  close(pb)
  
  out <- list(df = df, mosaics = list_mosaics, hotspots = hotspots, paths = paths, foraging_area = foraging_area)
  invisible(out)
}