# Gaura State-space
# Setup----
library(dplyr)
library(tidyr)
library(ggplot2)
library(rjags)
# library(runjags)
library(stringr)
library(dclone)
library(MCMCvis)
'%!in%' <- function(x,y)!('%in%'(x,y))
theme_set(theme_bw(base_size = 14))


# Read in input files
census <- read.csv("cobp.count.csv", header = TRUE)
reps <- read.csv("count.replication.csv")
# # by 13 creek segments (removing 2 with zeros in most years)
# ydata <- census %>%
#   # filter(segment %!in% c("C-VII", "C-VIII")) %>% 
#   spread(year, flower.count) %>% 
#   dplyr::select(-creek, -segment) %>%
#   data.matrix()
# 
# # by 3 creeks
# ydata <- census %>%
#   group_by(year, creek) %>%
#   summarise(flower.count = sum(flower.count)) %>% 
#   spread(year, flower.count) %>% 
#   dplyr::select(-creek) %>%
#   as.matrix()

segments <- read.csv("data/segment.counts.csv")
segdat <- segments %>% 
  dplyr::select(-X) %>% 
  gather(segment, flower.count, C.I:U.II) %>% 
  filter(is.na(year) == FALSE) %>% 
  mutate(segment = gsub(pattern = ".", replacement = "-", x = segment, fixed = TRUE))

ydata <- segdat %>%
  filter(segment %!in% c("C-VIII")) %>%
  spread(year, flower.count) %>% 
  dplyr::select(-segment) %>%
  data.matrix()

polygons <- read.csv("data/polygon.counts.csv")
polydat <- polygons %>% 
  gather(year, flower.count, X2002:X2018) %>% 
  mutate(year = gsub(pattern = "X", replacement = "", x = year, fixed = TRUE),
         flower.count = as.numeric(gsub(pattern = " ", replacement = "", x = as.character(flower.count), fixed = TRUE))) %>% 
  group_by(segment, year) %>% 
  summarise(flower.count = sum(flower.count, na.rm = TRUE)) %>% 
  ungroup() %>% 
  filter(segment %!in% c("C-VIII")) %>%
  spread(year, flower.count) %>% 
  dplyr::select(-segment) %>%
  data.matrix()

ydata[,14:30] <- polydat # some excel errors noted in segment sums compared to polygon totals

# make NA zero, even if not quite right
ydata[is.na(ydata)] <- 0

# covariates
temp <- read.csv("temp.cheyenne.csv", header = TRUE, skip = 1) %>% 
  dplyr::select(year, temp.summer.apr.aug) %>% 
  filter(year > 1988) %>% 
  mutate(ztemp = scale(temp.summer.apr.aug))
prec <- read.csv("precip.cheyenne.csv", header = TRUE, skip = 1) %>% 
  dplyr::select(year, precip.summer.apr.aug) %>% 
  filter(year > 1988) %>% 
  mutate(zprec = scale(precip.summer.apr.aug))


# Settings
year1 <- min(census$year) # First year of estimation data
year2 <- max(census$year) # Last year of estimation data/first year of projection
year3 <- 2020 # Last year of projection
# thresholds <- c(1,10,50) # Quasi-extirpation levels to monitor

nYears <- year3 - year1 + 1
nPast <- year2 - year1 + 1
nFuture <- year3 - year2

zyear <- as.numeric(scale(year1:year3))

nSites <- nrow(ydata)

# From Chandler & Hostettler
# For projection, fill in the future counts with all NA
ydata <- cbind(ydata, matrix(NA, nrow=nSites, ncol=nFuture))
ztemp <- c(temp$ztemp, rep(0, nFuture))
zprec <- c(prec$zprec, rep(0, nFuture))


fakeN <- round((ydata + 1) * 2) # Try to set up plausible initial values for N
fillNAs <- function(N, default) { # Need to replace NAs in initial N with values
  tMax <- ncol(N)
  defaults <- which(is.na(N[,1]) & is.na(N[,2]))
  N[defaults,1] <- default # Assign missing values from the 1st year which also have missing values for the 2nd year with default supplied
  subsequent <- which(is.na(N[,1]))
  N[subsequent,1] <- N[subsequent,2] # Assign missing values from the 1st year which don't have missing values for the 2nd year with the 2nd year's value
  for (t in 2:(tMax-1)) {
    # Assign missing values from this year which also have missing values for the next year with the values from the previous year
    previous <- which(is.na(N[,t]) & is.na(N[,t+1]))
    N[previous,t] <- N[previous,t-1]
    # Assign missing values from this year which don't have missing values for the next year with the mean of the previous and next year
    means <- which(is.na(N[,t]))
    N[means,t] <- ceiling((N[means, t-1] + N[means, t+1])/2)
  }
  # Assign missing values from the last year with the values from the previous year
  previous <- which(is.na(N[,tMax]))
  N[previous, tMax] <- N[previous, tMax - 1]
  return(N)
}
fakeN <- fillNAs(fakeN, 1)
head(fakeN)


# State space Chandler Hostettler try 2----

# ricker model, det prob varies by site within yearly variation
# annual env. variation in R common to all sites
# rmax varies by site, K varies by site

# K way too large in below model, not sure why, poorer fit than above
# probably variation in counts can be explained by too many annual vars: R and p

sink(file="gaura.dm.rosette2.txt")
cat("
    # Regional stochasticity
    # Negative binomial initial abundance, Ricker-logistic dynamics + immigration
    model {
    lambda ~ dunif(0, 5000)
    alpha ~ dunif(0, 20)
    P <- alpha/(alpha + lambda)
    iota ~ dunif(0, 10)
    
    # growth rate 
    mean.r ~ dnorm(0, 0.001)
    
    sigma.nu ~ dunif(0, 5)
    tau.nu <- 1 / (sigma.nu * sigma.nu)
    
    
    mean.flower ~ dnorm(-2, .5)T(-5,-1)
    sigma.flower  ~ dunif(0, 2) 
    tau.flower <- 1 / (sigma.flower * sigma.flower)    
    
    
    # varying slopes by site
    for (v in 1:10){
    beta.mean[v] ~ dunif(-5, 5)
    sigma.beta[v] ~ dunif(0, 2)
    tau.beta[v] <- 1 / (sigma.beta[v] * sigma.beta[v])
    for (i in 1:nSites){
    site.beta[i, v] ~ dnorm(beta.mean[v], tau.beta[v])
    }
    }
    
    
    
    # vital rates
    for (i in 1:nSites){
    logit(flower[i,1]) <- norm.flower[i,1]
    norm.flower[i,1] ~ dnorm(mean.flower, tau.flower)T(-5,-1)
    
    for(t in 2:nYears) {
    
      log(lam[i,t-1]) <- norm.r[i,t-1]
      norm.r[i,t-1] ~ dnorm(r[i,t-1], tau.nu)T(-3,3)
      r[i,t-1] <-     mean.r +      
                      site.beta[i,1] * zprec[t-1] + 
                      site.beta[i,2] * zprec[t] + 
                      site.beta[i,3] * ztemp[t-1] + 
                      site.beta[i,4] * ztemp[t] +
                      site.beta[i,5] * log(N[i,t-1] + 1)
    
      logit(flower[i,t]) <- norm.flower[i,t] 
      norm.flower[i,t] ~ dnorm(link.flower[i,t], tau.flower)T(-5,-1)
      link.flower[i,t] <- mean.flower  +  
                      site.beta[i,6] * zprec[t-1] + 
                      site.beta[i,7] * zprec[t] + 
                      site.beta[i,8] * ztemp[t-1] + 
                      site.beta[i,9] * ztemp[t] +
                      site.beta[i,10] * log(N[i,t-1] + 1)
    
    }
    } 
    
    
    # Process model
    
    for(i in 1:nSites) {
    N[i,1] ~ dnegbin(P, alpha)
    for(t in 2:nYears){
    N[i,t] ~ dpois(muN[i, t-1])
    muN[i,t-1] <- N[i, t-1] * lam[i,t-1] + iota
    }
    
    for(t in 1:nYears) {
    # y[i,t] ~ dbin(0.87, N.fl[i,t])
    # y[i,t] ~ dpois(N.fl[i,t])
    y[i,t] ~ dbin(flower[i,t], N[i,t])
    }
    }
    
    # derived predictions summed across creeks
    for (t in 1:nYears){
    creek[,t] <- c(sum(N[c(1:7), t]), sum(N[c(8:12), t]), sum(N[c(13:14), t]))
    total[t] <- sum(creek[,t])
    }
    
    
    }
    ", fill=TRUE)
sink()

dat.rosette <- list(nSites=nSites, nYears=nYears, y=ydata, ztemp = ztemp, zprec = zprec, zyear = zyear,
                    creek = matrix(data = NA, nrow = 3, ncol = nYears),
                    total = rep(NA, nYears))
init.rosette <- function() list(lambda=runif(1, 0, 5000), alpha=runif(1, 0, 20), 
                                r=rnorm(nSites, mean = 1, sd = 1), 
                                mean.r = rnorm(1), sigma.r = runif(1, 0, 5),
                                K=rpois(nSites, apply(fakeN, MARGIN = 1, max)), 
                                iota=runif(1, 0, 10),
                                N=fakeN, 
                                sigma.nu = runif(1, 0, 5),
                                fl.mean = rnorm(1, -1.5, .7), # strong prior on fl.mean 
                                sigma.fl = runif(1, 0, 1),
                                beta1 = rnorm(1), beta2 = rnorm(1))
# sigma.fl = runif(1, 0, 1))
pars.rosette <- c("lambda", "alpha", "r", "sigma.nu", "iota", "N", "mean.flower", "flower", "site.beta", "lam", 
                  "sigma.beta", "beta.mean", "mean.r", "creek", "total")

jm.rosette <- jags.model("gaura.dm.rosette2.txt", dat.rosette, n.chains=3, n.adapt=100000)
update(jm.rosette, n.iter=200000)
jc.rosette <- coda.samples(jm.rosette, pars.rosette, n.iter=10000, thin = 10)
saveRDS(jc.rosette, "rosette2_mcmc.rds")

#100000 (50000 adapt)
#1000 sample

# check MCMC ----
summary(jc.rosette)
windows(record=T)
plot(jc.rosette)
warnings()
gelman.diag(jc.rosette)

jc <- jc.rosette
# gelman diagnostics for each variable separately
g <- matrix(NA, nrow=nvar(jc), ncol=2)
for (v in 1:nvar(jc)) {
  g[v,] <- gelman.diag(jc[,v])$psrf
}

MCMCsummary(jc)
# Dotplot of covariate effects
MCMCplot(jc.rosette, params = 'beta.mean',ref_ovl = TRUE,
         labels = c(paste("Population growth:", c('Precip (t-1)', 'Precip (t)', 
                                                  'Temp (t-1)', 'Temp (t)', 'N (t-1)'), sep = " "),
                    paste("Flowering rate:", c('Precip (t-1)', 'Precip (t)', 
                                               'Temp (t-1)', 'Temp (t)', 'N (t-1)'), sep = " ")))



# plot uncertainty ----

# could also grab rows from mcmclist, sort by variables, and plot population trajectories over time

ch <- rbind(jc.rosette[[1]], jc.rosette[[2]], jc.rosette[[3]])
vars <- dimnames(ch)[[2]]
# med <- apply(ch, MARGIN = 2, median)
qs <- apply(ch, MARGIN = 2, quantile, probs = c(0.25, .5, 0.75)) 
df <- data.frame(t(qs))
df$vars <- vars

ns <- df[grep(df$vars, pattern = "N[", fixed = TRUE),] %>% 
  mutate(segment = as.numeric(gsub("\\D+", "", str_split_fixed(vars, ",", 2)[,1])),
         year = 1989 - 1 + as.numeric(gsub("\\D+", "", str_split_fixed(vars, ",", 2)[,2])))

segs <- census %>%
  filter(segment %!in% c("C-VIII")) %>% 
  droplevels()

segnames <- levels(segs$segment)
ns$segment <- segnames[ns$segment]

counts_cr <- segs %>% 
  group_by(creek, year) %>% 
  summarise(flower.count = sum(flower.count)) %>% 
  mutate(segment = creek)
counts_all <- counts_cr %>% 
  group_by(year) %>% 
  summarise(flower.count = sum(flower.count),
            segment = "all")
counts <- bind_rows(counts_all, counts_cr, segs)

creek <- df[grep(df$vars, pattern = "creek", fixed = TRUE),] %>% 
  mutate(segment = rep(c("crow", "diamond", "unnamed"), nYears),
         year = 1989 - 1 + as.numeric(gsub("\\D+", "", str_split_fixed(vars, ",", 2)[,2])))

total <- df[grep(df$vars, pattern = "total", fixed = TRUE),] %>% 
  mutate(segment = c("all"),
         year = 1989:year3 )

ns <- bind_rows(ns, creek, total)

rs <- df[grep(df$vars, pattern = "r[", fixed = TRUE),]

bs <- df[grep(df$vars, pattern = "beta", fixed = TRUE),]

rs <- df[grep(df$vars, pattern = "flower[", fixed = TRUE),]


plt <- ggplot(data = ns %>% filter(year < 2026), aes(x = year, y = X50.)) +
  geom_point() +
  geom_linerange(aes(ymin = X25., ymax = X75.), alpha = 0.5) +
  geom_point(data = counts %>% filter(year < 2026), aes(x = year, y = flower.count), color = "red") +
  facet_wrap(~segment, scales = "free", ncol = 3) +
  ylab("Count") +
  ggtitle("Median prediction of total population size (with 50% CI) versus flower census")
plt  




