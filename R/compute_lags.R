# lag inference (from Athos) ####

compute_lag <- function(.xs, .skip_start=0, .skip_end=1, .scaling_factor=1e-2, .plot=FALSE) {
  if (length(.xs) < 3+.skip_start+.skip_end ) return(list())
  if (.skip_end<1) .skip_end <- 1
  
  UnifLagInt <- function(.T, .tau, .z, .delta){
    # use a uniform over delta  and integrate over lambda using Laplace, 
    # take the max of the log likelyhood 
    delta0 <- (.z%*%.delta) / (.delta%*%.delta);
    fdelta0 <- ((.z-delta0*.delta) %*% (.z-delta0*.delta)) / .T;
    F2Der <- ((-2*.delta%*%.delta)/.T) / fdelta0;
    # return(log(sqrt(pi/((1-T)*F2Der))) + ((1-T)/2*log(fdelta0)))
    return(-log((1-.T)*F2Der)+((1-.T)/2*log(fdelta0)) )
  }
  
  .xsc <- .xs # scaled variable
  .T <- length(.xs);
  .taus <- (1+.skip_start):(.T-.skip_end) # remove very unlikely taus (for speedup)
  .xsc <- .xsc * .scaling_factor
  .zs <- .xsc-mean(.xsc) # compute z-transformed var
  .DELTA_vecs <- lapply(.taus, function(.tau) 
    ifelse((1:.T)<.tau, -(.T-.tau)*(1+.T-.tau)/(2*.T), (1:.T) -.tau-(.T-.tau)*(1+.T-.tau)/(2*.T)) )
  .maxLikLag <- sapply(seq_along(.taus), function(.i) UnifLagInt(.T, .taus[.i], .zs, .DELTA_vecs[[.i]]))
  
  # browser()
  if (any(is.na(.maxLikLag))) {
    warning('some values of maxLikLag are NaN: you might want to increase the scaling factor')
    return(list())
  }
  .tau <- which.max(.maxLikLag) + .skip_start
  # case where the lag is longer than the data
  if (.tau > length(.maxLikLag)-.skip_end) {
    # plot(1:.T,.xs, col=ifelse((1:.T)>=.tau, "red", "black"))
    # lines(1:.T, c(mean(.xs)*rep(1,.T)))
    return(list(tau_idx=Inf, x_lag=mean(.xs), slope_after=NA))
  }
  .delta0 <- (.zs %*% .DELTA_vecs[[.tau-.skip_start]]) / (.DELTA_vecs[[.tau-.skip_start]] %*% .DELTA_vecs[[.tau-.skip_start]])
  .delta0 <- .delta0 / .scaling_factor
  .x0 <- mean(.xs) - .delta0 * (.T-.tau) * (1+.T-.tau)/(2*.T)
  
  if (.plot) {
    graphics::plot(1:.T,.xs, col=ifelse((1:.T)>=.tau, "red", "black"))
    graphics::lines(1:.T, c(.x0*rep(1,.tau), .x0+seq(1,.T-.tau)*.delta0))
    # .pl <- qplot(1:.T,.xs, col=(1:.T)>=.tau) +
    #   geom_line(aes(y=c(.x0*rep(1,.tau), .x0+seq(1,.T-.tau)*.delta0), col=NULL)) 
    # print(.pl)
    # return(list(tau_idx=.tau, x_lag=.x0, slope_after=.delta0, plot=.pl))
  }
  
  return(list(tau_idx=.tau, x_lag=.x0, slope_after=.delta0))
}

