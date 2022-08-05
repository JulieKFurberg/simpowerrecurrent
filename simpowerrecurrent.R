########################################################################################
### R function for simulation of recurrent events in the presence of terminal events ###
########################################################################################

# Required packages
require(survival)
require(mets)
require(ggplot2)

# Make this run without the preprocessing. 
setwd("~/Private/ErhvervsPhd/Article 4 - Simulation based sample size")
source("test_from_mets.R")

#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#
# Simulation of a single data set 
#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#

# - n: number of subjects
# - beta: treatment effect on recurrent events (Ghosh and Lin 2002)
# - gamma: treatment effect on death (Cox 1972)
# - mu0: baseline marginal mean function, mu_0(t)
# - Lam0D: baseline cumulative hazard of death, Lambda0D(t)
# - crate: constant rate of censoring, Lambda^C(t \mid Z) = Lambda_0^C(t) = 0.01
# - accrualtime: uniform accrual from 0 to accrualtime (should match timing of mu_0(t), Lambda0D(t))
# - admincens: administrative censoring at time adminscens (should match timing of mu_0(t), Lambda0D(t))

# Comments: Binary treatment variable, Z \sim bin(0.5, 1)



simRecurEnrollAdmin <- function(n = 10000, 
                                beta = 0,
                                gamma = 0,
                                mu0, 
                                Lam0D, 
                                crate = 0.01, 
                                accrualtime = 0, 
                                admincens = 8000){
  
  # n = 1000
  # beta = beta
  # gamma = gamma
  # mu0 = cumhaz_mu
  # Lam0D = cumhaz_S
  # crate = coef(linear_c)
  # accrualtime = 300 # in days
  # admincens = 1600
  
  
  dat <- simMarginalMeanCox(n,  
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
  
  return(updated)
}

## Example
mu0 <- function(time, a){a  * time}
S0 <- function(time, b){b * time}


times <- seq(0, 100, by = 1)

cumhaz_mu <- data.frame(time = times, 
                        cumhaz = mu0(time = times, a = 0.1))

cumhaz_S <- data.frame(time = times, 
                       cumhaz = S0(time = times, b = 0.01))
# Visual 
ggplot(aes(x = time, y = cumhaz), data = cumhaz_mu) + geom_line() + ggtitle("Mu_0(t) for recurrent events")
ggplot(aes(x = time, y = cumhaz), data = cumhaz_S) + geom_line() + ggtitle("Lambda_0^D(t) for death")


sim1 <- simRecurEnrollAdmin(n = 1000, 
                            beta = 0,
                            gamma = 0, 
                            mu0 = cumhaz_mu, 
                            Lam0D = cumhaz_S, 
                            crate = 0.01, 
                            accrualtime = 0, 
                            admincens = 8000)

head(sim1)


#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#
# Repeated simulation for power calculations
#--------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------#

# nsims = number of simulations
# alpha = significance level (0.05 default)

montecarlofunction <- function(nsims = 1000, 
                                n = 1000, 
                                beta =  0,
                                gamma = 0,
                                mu0, 
                                Lam0D, 
                                crate = 0.01, 
                                accrualtime = 0, 
                                admincens = 8000, 
                                alpha = 0.05){
  
  
  # n = 1000
  # beta = beta
  # gamma = gamma
  # mu0 = cumhaz_mu
  # Lam0D = cumhaz_S
  # crate = coef(linear_c)
  # accrualtime = 300 # in days
  # admincens = 1600
  
  r <- list()
  regmod <- list()
  resmat <- matrix(NA, nrow = nsims, ncol = 4) 
  colnames(resmat) <- c("beta", "sebeta", "reject?", "pval")
  
  for (i in 1:nsims){
    r[[i]] <- simRecurEnrollAdmin(n = n, 
                                  beta = beta,
                                  gamma = gamma,
                                  mu0 = mu0, 
                                  Lam0D = Lam0D, 
                                  crate = crate, 
                                  accrualtime = accrualtime, 
                                  admincens = admincens)
    # Overview of single data set
    dtable(r[[i]], ~ statusG + status + death, level = 2, response = 1)
    
    # Result from fitting GL regression model 
    regmod[[i]] <- recreg(Event(start, stop, statusG) ~ X1 + cluster(id), 
                          data = r[[i]], cause = 1, death.code = 3, cens.code = 0)
    
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


# Simulation 1000 trials
simres1 <- montecarlofunction(nsims = 100, n = 100, mu0 = cumhaz_mu, Lam0D = cumhaz_S, alpha = 0.05/2)

simres1

