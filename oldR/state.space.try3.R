#Gaura state-space models
#4/9/14

#updated 4/25 to try and incorporate all 13 subunits of subpopulations
#trying out binomial mix approach instead


#Three unconnected subpopulations of a rare plant
#Flowers counted every year, assumed to be metric for total population
#However, proportion of population in flowering stage varies.
#Also, annual variation in growth rate/proportion flowering correlated between subpopulations

library(rjags)
library(R2jags)

setwd("C:/Users/Tyson/Dropbox/Gaura")

counts <- read.csv("subcounts.csv", header = TRUE)
counts <- matrix(counts[,4], nrow = 25) #treating all subunits as unique populations, ignoring creek (should change later)
#Columns 1-2: Unnamed, 3:7 Diamond, 8:13 Crow 
counts <- counts[,c(1:13)]

veg <- read.csv("veg.csv", header = TRUE)
veg$ratio <- veg$Flow04/(veg$Flow04 + veg$Veg04)
veg$subdiv <- paste(veg$Creek, veg$Subunit, sep = ".")


climate <- read.csv("climate.csv", header = TRUE)
precip <- climate[, c("DATE", "TPCP")]
precip$DATE <- ymd(precip$DATE)
precip <- precip[year(precip$DATE) > 1985, ]
precip <- precip[month(precip$DATE) > 3, ]
precip <- precip[month(precip$DATE) < 9, ]

rain <- ddply(precip, .(year(DATE)), summarize,
              year = year(DATE)[1],
              mm = mean(TPCP)/10)

rain <- rain[, c(2:3)]
############################################
#State-space with varying proportion of population flowering


# Specify model in BUGS language
sink("model.txt")
cat("
    model {
    # Priors and constraints
    for (j in 1:13){
    alpha.lam[j] ~ dnorm(7, 0.01)T(1, 15)         #each population has different size
    trend[j] ~ dnorm(0, .01)
    }

#     a0 ~ dnorm(0, 0.01)
#     a1 ~ dnorm(0, 0.01)
#     a2 ~ dnorm(0, 0.01)
# 
#     b0 ~ dnorm(0, 0.01)
#     b1 ~ dnorm(0, 0.01)

#     error.a ~ dunif(0,10)
#     error.b ~ dunif(0,10)
  
    
    # Likelihood
    # State process
    for (t in 1:nyear){
      for (j in 1:13){
        N.est[t,j] ~ dpois(lambda[t,j])
        log(lambda[t,j]) <- alpha.lam[j] + trend[j] * t
#                             + a0*rain0[t,1] + a1*rain1[t,1] + a2*rain2[t,1]
#                             + b0*log(X1[t,j]) + b1*log(X2[t,j])
      }
    }
    
    # Observation process
    for (t in 1:nyear) {
      for (j in 1:13){
        y[t,j] ~ dbin(error[t,j], N.est[t,j])
        error[t,j] ~ dbeta(3,5)
      }
    }
    
    }
    ",fill = TRUE)
sink()

pyears <- 0 # Number of future years with predictions
y <- rbind(counts, matrix(data = rep(NA, pyears*13), ncol = 13))
y[3,7] <- 20 #replace NA with reasonable number
y[19,12] <- 10
year <- 1991:(2013 + pyears)
nyear  <- length(year)

# Bundle data
bugs.data <- list(y = y[c(3:25), ] + 1, nyear = nyear) 
#                   rain0 = scale(rain$mm[c(6:28)]), 
#                   rain1 = scale(rain$mm[c(5:27)]), 
#                   rain2 = scale(rain$mm[c(4:26)]),
#                   X1 = y[c(2:24), ], X2 = y[c(1:23), ])

# Initial values
inits <- function(){list(alpha.lam = log(5*y[1,]),
#                          error.a = runif(1, 1, 5), error.b = runif(1, 1, 5), 
                         trend = rnorm(13))}
#                          a0 = rnorm(1, 0, 1), a1 = rnorm(1, 0, 1), a2 = rnorm(1, 0, 1),
#                           b0 = rnorm(1, 0, 1), b1 = rnorm(1, 0, 1))}


# Parameters monitored
params <- c("alpha.lam", "error.a", "error.b", "trend") 

# MCMC settings
ni <- 5000
nt <- 100
nb <- 1000
nc <- 3

# Call WinBUGS from R (BRT 4 min)
out <- jags(bugs.data, inits, params, "model.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, digits = 2)


# Draw figure
fitted <- lower <- upper <- matrix(data = rep(NA, nyear*13), nrow = nyear, ncol = 13)
for (i in 1:nyear){
  for (j in 1:13){
    fitted[i,j] <- exp(mean(out$BUGSoutput$sims.list$logN.est[,i,j]))
    lower[i,j] <- exp(quantile(out$BUGSoutput$sims.list$logN.est[,i,j], 0.025))
    upper[i,j] <- exp(quantile(out$BUGSoutput$sims.list$logN.est[,i,j], 0.975))
  }
}

for (j in 1:13){
  m1 <- min(c(fitted[,j], y[,j], lower[,j]), na.rm = TRUE)
  m2 <- max(c(fitted[,j], y[,j], upper[,j]), na.rm = TRUE)
  par(mar = c(4.5, 4, 1, 1))
  plot(0, 0, ylim = c(m1, m2), xlim = c(1, nyear), ylab = "Population size", xlab = "Year", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  axis(2, las = 1)
  axis(1, at = 1:nyear, labels = year)
  polygon(x = c(1:nyear, nyear:1), y = c(lower[,j], upper[nyear:1,j]), col = "gray90", border = "gray90")
  points(y[,j], type = "l", col = "black", lwd = 2)
  points(fitted[,j], type = "l", col = "blue", lwd = 2)
}


# Probability of N(2015) < N(2009)
mean(out$BUGSoutput$sims.list$N.est[,26] < out$BUGSoutput$mean$N.est[20])