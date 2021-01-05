#Gaura state-space models
#4/9/14

#Three unconnected subpopulations of a rare plant
#Flowers counted every year, assumed to be metric for total population
#However, proportion of population in flowering stage varies.
#Also, annual variation in growth rate/proportion flowering correlated between subpopulations

library(rjags)
library(R2jags)

setwd("C:/Users/Tyson/Dropbox/Gaura")

counts <- read.csv("gaura.count.csv", header = TRUE)
counts <- counts[-c(1:2), ]


############################################
#State-space with varying proportion of population flowering


# Specify model in BUGS language
sink("model.txt")
cat("
    model {
# Priors and constraints

    N.est[1] ~ dnorm(1000, 0.0001)T(0,)       # Prior for initial population size
    K ~ dnorm(10000, 0.0001)T(0,)             # Prior for carrying capacity
  #mean demographic parameters
    l.p ~ dunif(0, 0.5)                     # Here, detection probability is proportion of plants flowering that year
    l.r ~ dnorm(0, 0.001) T(-5,5)     # Prior for mean growth rate

  #precision of standard deviations of temporal variability
    sigma.p ~ dunif(0, 1)
    tau.p <- pow(sigma.p, -2)

    sigma.r ~ dunif(0, 1)            # Prior for sd of mean growth rate
    tau.r <- pow(sigma.r, -2)
  
  #distribution of error terms (bounded to help with convergence)
    for (t in 1:(T-1)){
      eps.r[t] ~ dnorm(0, tau.r)T(-5,5)
    }

    for (t in 1:T){
      eps.p[t] ~ dnorm(0, tau.r)T(-5,5)
    }

#Constrain parameters
  for (t in 1:(T-1)){
    log(r[t]) <- l.r + eps.r[t]
  }
  for (t in 1:T){
    logit(p[t]) <- l.p + eps.p[t]
  }

#derived parameters
  mean.p <- exp(l.p)/(1+exp(l.p))
  mean.r <- exp(l.r)

  #population growth rate
  for (t in 1:(T-1)){
    lambda[t] <- N.est[t+1]/N.est[t]
  }

    # Likelihood
    # State process
    for (t in 1:(T-1)){
    N.est[t+1] <- N.est[t]*(1 + dens.dep.r[t])
    dens.dep.r[t] <- r[t]*(1 - (N.est[t] / K))
    }

    # Observation process
    for (t in 1:T) {
    y[t] ~ dbin(p[t], N.dis[t])
    N.dis[t] <- round(N.est[t])
    }
    
    }
    ",fill = TRUE)
sink()

# House martin population data from Magden
pyears <- 0 # Number of future years with predictions
hm <- c(counts$diamond, rep(NA, pyears))
year <- 1988:(2012 + pyears)

# Bundle data
bugs.data <- list(y = hm, T = length(year))

# Initial values
inits <- function(){list(N.est = c(runif(1, 2000, 12000), rep(NA, (length(year)-1))), K = runif(1, 1000, 100000),
                         l.p = runif(1, 0, 0.5), l.r = runif(1, -3, 3),
                         sigma.p = runif(1, 0, 1), sigma.r = runif (1, 0, 1))}

# Parameters monitored
params <- c("N.est", "K", "mean.p", "mean.r", "p", "r", "sigma.p", "sigma.r", "lambda") 

# MCMC settings
ni <- 1000
nt <- 10
nb <- 200
nc <- 3

# Call WinBUGS from R (BRT 4 min)
out <- jags(bugs.data, inits, params, "model.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, digits = 1)


# Draw figure
fitted <- lower <- upper <- numeric()
n.years <- length(hm)
for (i in 1:n.years){
  fitted[i] <- mean(out$BUGSoutput$sims.list$N.est[,i])
  lower[i] <- quantile(out$BUGSoutput$sims.list$N.est[,i], 0.025)
  upper[i] <- quantile(out$BUGSoutput$sims.list$N.est[,i], 0.975)}
m1 <- min(c(fitted, hm, lower), na.rm = TRUE)
m2 <- max(c(fitted, hm, upper), na.rm = TRUE)
par(mar = c(4.5, 4, 1, 1))
plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Population size", xlab = "Year", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
axis(2, las = 1)
axis(1, at = 1:n.years, labels = year)
polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
points(hm, type = "l", col = "black", lwd = 2)
points(fitted, type = "l", col = "blue", lwd = 2)
legend(x = 1, y = 150, legend = c("Counts", "Estimates"), lty = c(1, 1), lwd = c(2, 2), col = c("black", "blue"), bty = "n", cex = 1)

# Probability of N(2015) < N(2009)
mean(out$BUGSoutput$sims.list$N.est[,26] < out$BUGSoutput$mean$N.est[20])










#######################################
#Simplest state-space, density independent, one plant population

# Specify model in BUGS language
sink("model.txt")
cat("
    model {
    # Priors and constraints
    logN.est[1] ~ dnorm(8, 0.01)       # Prior for initial population size
    mean.r ~ dnorm(1, 0.001)             # Prior for mean growth rate
    sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
    sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    r[t] ~ dnorm(mean.r, tau.proc)
    logN.est[t+1] <- logN.est[t] + r[t]
    }
    # Observation process
    for (t in 1:T) {
    y[t] ~ dnorm(logN.est[t], tau.obs)
    }
    
    # Population sizes on real scale
    for (t in 1:T) {
    N.est[t] <- exp(logN.est[t])
    }
    }
    ",fill = TRUE)
sink()

# House martin population data from Magden
pyears <- 10 # Number of future years with predictions
hm <- c(counts$diamond, rep(NA, pyears))
year <- 1986:(2012 + pyears)

# Bundle data
bugs.data <- list(y = log(hm), T = length(year))

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 1), mean.r = rnorm(1), sigma.obs = runif(1, 0, 1), logN.est = c(rnorm(1, 8, 0.1), rep(NA, (length(year)-1))))}

# Parameters monitored
params <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 50000
nt <- 6
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT 3 min)
#ssm <- bugs(bugs.data, inits, parameters, "ssm.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Call WinBUGS from R (BRT 4 min)
out <- jags(bugs.data, inits, params, "model.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, digits = 1)


# Draw figure
fitted <- lower <- upper <- numeric()
n.years <- length(hm)
for (i in 1:n.years){
  fitted[i] <- mean(out$BUGSoutput$sims.list$N.est[,i])
  lower[i] <- quantile(out$BUGSoutput$sims.list$N.est[,i], 0.025)
  upper[i] <- quantile(out$BUGSoutput$sims.list$N.est[,i], 0.975)}
m1 <- min(c(fitted, hm, lower), na.rm = TRUE)
m2 <- max(c(fitted, hm, upper), na.rm = TRUE)
par(mar = c(4.5, 4, 1, 1))
plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Population size", xlab = "Year", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
axis(2, las = 1)
axis(1, at = 1:n.years, labels = year)
polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
points(hm, type = "l", col = "black", lwd = 2)
points(fitted, type = "l", col = "blue", lwd = 2)
legend(x = 1, y = 150, legend = c("Counts", "Estimates"), lty = c(1, 1), lwd = c(2, 2), col = c("black", "blue"), bty = "n", cex = 1)

# Probability of N(2015) < N(2009)
mean(out$BUGSoutput$sims.list$N.est[,26] < out$BUGSoutput$mean$N.est[20])



########################################################################
#
# 5. State-space models
#
###########################################################################

# 5.2. A simple model
n.years <- 25           # Number of years
N1 <- 30                # Initial population size
mean.lambda <- 1.02     # Mean annual population growth rate
sigma2.lambda <- 0.02   # Process (temporal) variation of the growth rate
sigma2.y <- 20          # Variance of the observation error

y <- N <- numeric(n.years)
N[1] <- N1
lambda <- rnorm(n.years-1, mean.lambda, sqrt(sigma2.lambda))
for (t in 1:(n.years-1)){
  N[t+1] <- N[t] * lambda[t]
}

for (t in 1:n.years){
  y[t] <- rnorm(1, N[t], sqrt(sigma2.y))
}

# Specify model in BUGS language
sink("ssm.bug")
cat("
    model { 
    # Priors and constraints
    N.est[1] ~ dunif(0, 500)            # Prior for initial population size
    mean.lambda ~ dunif(0, 10)          # Prior for mean growth rate
    sigma.proc ~ dunif(0, 10)           # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    sigma.obs ~ dunif(0, 100)           # Prior for sd of observation process
    sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    lambda[t] ~ dnorm(mean.lambda, tau.proc) 
    N.est[t+1] <- N.est[t] * lambda[t] 
    }
    # Observation process
    for (t in 1:T) {
    y[t] ~ dnorm(N.est[t], tau.obs)
    }
    }
    ",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(y = y, T = n.years)

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 5), mean.lambda = runif(1, 0.1, 2), sigma.obs = runif(1, 0, 10), N.est = c(runif(1, 20, 40), rep(NA, (n.years-1))))} 

# Parameters monitored
parameters <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT <1 min)
ssm <- bugs(bugs.data, inits, parameters, "ssm.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Define function to draw a graph to summarize results
graph.ssm <- function(ssm, N, y){
  fitted <- lower <- upper <- numeric()
  n.years <- length(y)
  for (i in 1:n.years){
    fitted[i] <- mean(ssm$sims.list$N.est[,i])
    lower[i] <- quantile(ssm$sims.list$N.est[,i], 0.025)
    upper[i] <- quantile(ssm$sims.list$N.est[,i], 0.975)}
  m1 <- min(c(y, fitted, N, lower))
  m2 <- max(c(y, fitted, N, upper))
  par(mar = c(4.5, 4, 1, 1), cex = 1.2)
  plot(0, 0, ylim = c(m1, m2), xlim = c(0.5, n.years), ylab = "Population size", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2, frame = FALSE, axes = FALSE)
  axis(2, las = 1)
  axis(1, at = seq(0, n.years, 5), labels = seq(0, n.years, 5))
  axis(1, at = 0:n.years, labels = rep("", n.years + 1), tcl = -0.25)
  polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
  points(N, type = "l", col = "red", lwd = 2)
  points(y, type = "l", col = "black", lwd = 2)
  points(fitted, type = "l", col = "blue", lwd = 2)
  legend(x = 1, y = m2, legend = c("True", "Observed", "Estimated"), lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("red", "black", "blue"), bty = "n", cex = 1)
}

# Execute function: Produce figure 
graph.ssm(ssm, N, y)


# 5.3. Systematic bias in the observation process
n.years <- 25  # Number of years
N <- rep(50, n.years) 

p <- 0.7
y <- numeric(n.years)
for (t in 1:n.years){
  y[t] <- rbinom(1, N[t], p)
}
y

# Bundle data
bugs.data <- list(y = y, T = n.years)

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 5), mean.lambda = runif(1, 0.1, 2), sigma.obs = runif(1, 0, 10), N.est = c(runif(1, 30, 60), rep(NA, (n.years-1))))}

# Parameters monitored
parameters <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 25000
nt <- 3
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT <1 min)
ssm <- bugs(bugs.data, inits, parameters, "ssm.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(ssm, digits = 3)

# Produce figure
graph.ssm(ssm, N, y)

n.years <- 25  # Number of years
N <- rep(50, n.years)

lp <- -0.5 + 0.1*(1:n.years)  # Increasing trend of logit p
p <- plogis(lp)
y <- numeric(n.years)
for (t in 1:n.years){
  y[t] <- rbinom(1, N[t], p[t])
}

# Bundle data
bugs.data <- list(y = y, T = n.years)

# Call WinBUGS from R (BRT <1 min)
ssm <- bugs(bugs.data, inits, parameters, "ssm.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Produce figure
graph.ssm(ssm, N, y)
points(N*p, col = "black", type = "l", lwd = 2, lty = 2)
legend(x = 1, y = 45.5, legend = "Np", lwd = 2, col = "black", lty = 2, bty = "n")


# 5.4. Real example: House martin population counts in the village of Magden
# Specify model in BUGS language
sink("ssm.bug")
cat("
    model {
    # Priors and constraints
    logN.est[1] ~ dnorm(5.6, 0.01)       # Prior for initial population size
    mean.r ~ dnorm(1, 0.001)             # Prior for mean growth rate
    sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
    sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    r[t] ~ dnorm(mean.r, tau.proc)
    logN.est[t+1] <- logN.est[t] + r[t]
    }
    # Observation process
    for (t in 1:T) {
    y[t] ~ dnorm(logN.est[t], tau.obs)
    }
    
    # Population sizes on real scale
    for (t in 1:T) {
    N.est[t] <- exp(logN.est[t])
    }
    }
    ",fill = TRUE)
sink()

# House martin population data from Magden
pyears <- 6 # Number of future years with predictions
hm <- c(271, 261, 309, 318, 231, 216, 208, 226, 195, 226, 233, 209, 226, 192, 191, 225, 245, 205, 191, 174, rep(NA, pyears))
year <- 1990:(2009 + pyears)

# Bundle data
bugs.data <- list(y = log(hm), T = length(year))

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 1), mean.r = rnorm(1), sigma.obs = runif(1, 0, 1), logN.est = c(rnorm(1, 5.6, 0.1), rep(NA, (length(year)-1))))}

# Parameters monitored
parameters <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 200000
nt <- 6
nb <- 100000
nc <- 3

# Call WinBUGS from R (BRT 3 min)
hm.ssm <- bugs(bugs.data, inits, parameters, "ssm.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(hm.ssm, digits = 3)

# Draw figure
fitted <- lower <- upper <- numeric()
year <- 1990:2015
n.years <- length(hm)
for (i in 1:n.years){
  fitted[i] <- mean(hm.ssm$sims.list$N.est[,i])
  lower[i] <- quantile(hm.ssm$sims.list$N.est[,i], 0.025)
  upper[i] <- quantile(hm.ssm$sims.list$N.est[,i], 0.975)}
m1 <- min(c(fitted, hm, lower), na.rm = TRUE)
m2 <- max(c(fitted, hm, upper), na.rm = TRUE)
par(mar = c(4.5, 4, 1, 1))
plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Population size", xlab = "Year", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
axis(2, las = 1)
axis(1, at = 1:n.years, labels = year)
polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
points(hm, type = "l", col = "black", lwd = 2)
points(fitted, type = "l", col = "blue", lwd = 2)
legend(x = 1, y = 150, legend = c("Counts", "Estimates"), lty = c(1, 1), lwd = c(2, 2), col = c("black", "blue"), bty = "n", cex = 1)

# Probability of N(2015) < N(2009)
mean(hm.ssm$sims.list$N.est[,26] < hm.ssm$mean$N.est[20])


