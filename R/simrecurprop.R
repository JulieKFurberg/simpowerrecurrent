#' Function that simulates a single data set according to the proportional means
#' and proportional hazards model
#' 
#' @param n Total sample size. Randomisation 1:1
#' @param beta Treatment effect on recurrent events
#' @param gamma Treatment effect on terminal events
#' @param mu0 Reference rate for marginal mean. Data frame with times and mu0(times)
#' @param Lam0D Reference rate for cumulative hazard of death. Data frame with times and Lam0D(times)
#' @param crate Exponential censoring rate during trial
#' @param accrualtime Max accrual time. Same time unit as mu0 and Lam0D
#' @param admincens Study closure, time of administrative censoring.  Same time unit as mu0 and Lam0D
#' 
#' @keywords simpowerrecurrent

#' @export
simrecurprop <- function(n,
                         beta = 0,
                         gamma = 0,
                         mu0,
                         Lam0D,
                         crate = 0.01,
                         accrualtime = 0,
                         admincens = 8000){
  
  dat <- mets:::simMarginalMeanCox(n,
                                   cens = crate,
                                   k1 = 1,
                                   Lam1 = mu0,
                                   LamD = Lam0D,
                                   beta1 = c(beta, 0),
                                   betad = c(gamma, 0))
  
  # Overview
  dtable(dat, ~ statusG + status + death, level = 2, response = 1)
  
  
  ###  Add uniform enrollment times  ###
  enroll <- runif(n, 0, accrualtime)
  
  # Move all start and stop dates per individual
  unisub <- length(unique(dat$id))
  idi <- list()
  for (i in 1:unisub){
    idi[[i]] <- subset(dat, id == unique(dat$id)[[i]])
    idi[[i]]$start <- idi[[i]]$start + enroll[i]
    idi[[i]]$stop <- idi[[i]]$stop + enroll[i]
  }
  
  dat_with_enroll <- do.call("rbind", idi)
  dat_with_enroll
  
  ### Make administrative censoring at time 'admincens' xxx
  
  if (any(dat_with_enroll$stop > admincens)){
    admincensdat <- subset(dat_with_enroll, stop > admincens)
    
    
    # First records post admin cens - new cens
    firstafteradmin <- admincensdat[!duplicated(admincensdat$id), ]
    
    firstafteradmin$stop <- admincens
    firstafteradmin$statusG <- 0
    
    # Collect
    notadmincens <- subset(dat_with_enroll, stop <= admincens)
    both <- rbind(notadmincens, firstafteradmin)
    sortboth <- both[order(both$id),]
    
  }
  else {sortboth <- dat_with_enroll}
  updated <- sortboth
  
  # Drop irrelevant variables
  updated$status <- ifelse(updated$statusG == 3, 2, updated$statusG)
  updated$Z <- updated$X1  
  updated <- updated[, c("id", "start", "stop", "status", "Z")]
  
  return(updated)
}