#' Function that makes power approximation based on simulated data using \code{simrecurprop}
#'
#' @param nsims Total number of simulations
#' @param n Total sample size. Randomisation 1:1
#' @param beta Treatment effect on recurrent events
#' @param gamma Treatment effect on terminal events
#' @param mu0 Reference rate for marginal mean. Data frame with times and mu0(times)
#' @param Lam0D Reference rate for cumulative hazard of death. Data frame with times and Lam0D(times)
#' @param crate Exponential censoring rate during trial
#' @param accrualtime Max accrual time. Same time unit as mu0 and Lam0D
#' @param admincens Study closure, time of administrative censoring.  Same time unit as mu0 and Lam0D
#' @param alpha Significance level. Default is 0.05
#' 
#' @keywords simpowerrecurrent

#' @export
#' 
powerest <- function(nsims = 1000,
                     n = 1000,
                     beta =  0,
                     gamma = 0,
                     mu0,
                     Lam0D,
                     crate = 0.01,
                     accrualtime = 0,
                     admincens = 8000,
                     alpha = 0.05){
  
  r <- list()
  regmod <- list()
  resmat <- matrix(NA, nrow = nsims, ncol = 4)
  colnames(resmat) <- c("beta", "sebeta", "reject?", "pval")
  
  for (i in 1:nsims){
    r[[i]] <- simrecurprop(n = n,
                           beta = beta,
                           gamma = gamma,
                           mu0 = mu0,
                           Lam0D = Lam0D,
                           crate = crate,
                           accrualtime = accrualtime,
                           admincens = admincens)
    
    # Result from fitting GL regression model
    regmod[[i]] <- recreg(Event(start, stop, status, cens = (status == 0)) ~ Z + cluster(id),
                          data = r[[i]], cause = 1, death.code = 2, cens.code = 0)
    
    # Collecting results of interest
    resmat[i,] <- c(regmod[[i]]$coef,
                    regmod[[i]]$se.coef,
                    as.numeric(abs(regmod[[i]]$coef / regmod[[i]]$se.coef) >= abs(qnorm(alpha/2))),
                    2*(1 - pnorm(abs(regmod[[i]]$coef / regmod[[i]]$se.coef)))
    )
    
  }
  
  
  # output
  out <- list(resmat = resmat,
              power = mean(resmat[,3]),
              betamean = mean(na.omit(resmat[,1])),
              betasemean = mean(na.omit(resmat[,2])),
              na.obs = sum(is.na(resmat[,3]))
  )
  return(out)
}