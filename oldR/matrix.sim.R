#simple matrix model simulation for JAGS analysis
setwd("C:/Users/Tyson/Dropbox/Gaura")


#gaura matrix model using floyd & ranker raw data
library(rjags)
library(R2jags)
library(popbio)
library(plyr)
library(reshape)


#plant with simpler dynamics
germ <- 0.2
surv.sd <- 0.7
surv.sm <- 0.7
psi.sm <- 0.3
psi.fl <- 1 - psi.sm
fert <- 10

n <- c(100, 50, 10, 0)
n.states <- 4

Trans <- matrix(c(
  surv.sd * (1 - germ), 0,                0, 0,
  surv.sd * germ,       surv.sm * psi.sm, 0, 0,
  0,                    surv.sm * psi.fl, 0, 0,
  1 - surv.sd,          1 - surv.sm,      1, 1), nrow = n.states, byrow = TRUE)

#use multinomial draws to project population structure in next census
proj.n <- function(n, Trans, fert, surv.sd){
  n1 <- matrix(NA, nrow = 4, ncol = 4)
  seeds <- n[1] + n[3] * rpois(1, fert) 
  n1[,1] <- rmultinom(1, seeds, prob = Trans[,1])  
  n1[,2] <- rmultinom(1, n[2], prob = Trans[,2])  
  n1[,3] <- rmultinom(1, n[3], prob = Trans[,3])
  n1[,4] <- rmultinom(1, n[4], prob = Trans[,4])
  return(n1)
}

#simulate 20 transitions, then see if JAGS can figure out the parameters

pop.1 <- pop.2 <- pop.3 <- matrix(NA, nrow = 4, ncol = 20)
trans.1  <- trans.2 <- array(data = NA, dim = c(4, 4, 20))

for (i in 1:20){
  pop.1[,i] <- c(rpois(1, 100), rpois(1, 50), rpois(1,10), 0)
  trans.1[,,i] <- proj.n(pop.1[,i], Trans, fert, surv.sd)
  pop.2[,i] <- rowSums(trans.1[,,i])
  pop.2[4,i] <- 0
  trans.2[,,i] <- proj.n(pop.2[,i], Trans, fert, surv.sd)
  pop.3[,i] <- rowSums(trans.2[,,i])
  pop.3[4,i] <- 0
}   



#specify model in JAGS language
sink("matrix.model.txt")
cat("
    model{
    
    #Parameters
    #surv: Survival for each stage
    #psi: transitions from each stage
    #fert: fertility/flower
    #germ: recruitment rate from seed
    
    #States:
    #     1 seed
    #     2 small
    #     3 flower
    #     4 dead
    
    #Priors and constraints
    #Survival for each stage
    surv.sd ~ dunif(0,1)
    surv.sm ~ dunif(0,1)
    
    #Transitions
    psi.sm ~ dunif(0,1)
    psi.fl <- 1 - psi.sm
    
    #Other parameters
    fert ~ dpois(10, .001)
    germ ~ dunif(0,1)

    #derived
    mean.sb <- mean(pop[1,])
    
    #Likelihood

    for (i in 1:40){
      new.seeds[i] ~ dpois(lambda[i])
      lambda[i] <- fert * pop[3,i]
      all.seeds[i] <- pop[1,i] + new.seeds[i]
      fate[1,1,i] <- all.seeds[i] * surv.sd * (1 - germ)
      fate[2,1,i] <- all.seeds[i] * surv.sd * germ
    
      fate[2,2,i] <- pop[2,i] * surv.sm * psi.sm
      fate[3,2,i] <- pop[2,i] * surv.sm * psi.fl
    }
    
    
    }", fill = TRUE)
sink()

#Bundle data

pop <- cbind(pop.1, pop.2)
pop <- pop[c(1:3),]
fate <- array(data = NA, dim = c(4, 4, 40))
fate[,,c(1:20)] <- trans.1
fate[,,c(21:40)] <- trans.2
fate <- fate[,c(1:3),]

bugs.data <- list(pop = pop, fate = fate)

#Initial values
inits <- function(){list(surv.sd = runif(1,0,1), 
                         surv.sm = runif(1,0,1),
                          psi.sm = runif(1,0,1),
                          fert = rnorm(1,10,1),
                          germ = runif(1,0,1))}

#Parameters monitored
params <- c("surv.sd", "surv.sm", "psi.sm", "psi.fl", "fert", "germ", "mean.sb")

#MCMC settings
ni <- 5000
nt <- 6
nb <- 2000
nc <- 3

#Call JAGS
out <- jags(bugs.data, inits, params, "matrix.model.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, digits = 2)


