# Gaura State-space
# Setup----
library(dplyr)
library(tidyr)
library(ggplot2)
library(rjags)
library(MCMCvis)
# library(runjags)
library(stringr)
library(dclone)
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


nSites <- nrow(ydata)
# nYears <- ncol(ydata)

# From Chandler & Hostettler
# For projection, fill in the future counts with all NA
ydata <- cbind(ydata, matrix(NA, nrow=nSites, ncol=nFuture))
ztemp <- c(temp$ztemp, rep(0, nFuture))
zprec <- c(prec$zprec, rep(0, nFuture))
# temporal trend
zyear <- as.numeric(scale(year1:year3))

fakeN <- round((ydata + 1) * 5) # Try to set up plausible initial values for N
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


# IPM----

# Two stage demographic matrices from Floyd & Ranker study
# Use as priors from Daniel's collapsed transition matrices from F&R data
# Will likely expand range of possible values since these are only from 2 years' transitions
# A[1,1] <- mean .7 (range .6 to .75)
# A[2,1] <- mean .137 (range .1 to .22)
# A[1,2] <- mean 9.7 (range 2.3 to 22.7) 
# A[2,2] <- mean 0.016 (range 0 to 0.05) (possibly ignore flower survival)

# Transition matrix components of survival, flowering, and germination
# Any of these components could be modified by temp/precip or density dependence
# A[1,1] <- survR * (1 - flower)
# A[2,1] <- survR * flower
# A[1,2] <-  germ * survR  # unknown really when seeds germinate or what role flower fecundity in t-1 plays
# A[2,2] <- survF # could make zero for simplicity


# i = site/segment
# t = year

# Basic demography (no covariates)
# Track rosettes and flowers as R[i, t] and F[i, t]
# R[i, t + 1] <- R[i,t] * A[1,1] + F[i,t] * A[1,2]
# F[i, t + 1] <- R[i,t] * A[2,1] + F[i,t] * A[2,2]

sink(file="gaura.ipm1.txt")
cat("
    model {
    # Negative binomial initial abundance
    lambda ~ dunif(0, 5000)
    alpha ~ dunif(0, 20)
    P <- alpha/(alpha + lambda)
    
    # Vital rates (annual variation)
    # mean and standard deviations
    # on link scales (log and logit)
    mean.survR ~ dnorm(1.5, .25)  # F & R rosette survival is 80 - 87% (including ones that then flower)
    mean.flower ~ dnorm(-2, .5) # based on rosette:flower ratios and flowering rates from F&R, 
    # survF[i] ~ dunif(0, 0.05)
    mean.germ ~ dunif(-3, 3.5) # combination of recruitment and survival from F&R (0 to 20 per flower)
    mean.iota ~ dpois(1) # immigration

    
    # sigma.site.survR ~ dunif(0, 5)
    # tau.site.survR <- 1 / (sigma.site.survR * sigma.site.survR)
    # sigma.site.flower ~ dunif(0, 5)
    # tau.site.flower <- 1 / (sigma.site.flower * sigma.site.flower)
    # sigma.site.germ ~ dunif(0, 5)
    # tau.site.germ <- 1 / (sigma.site.germ * sigma.site.germ)

    sigma.year.survR ~ dunif(0, 5)
    tau.year.survR <- 1 / (sigma.year.survR * sigma.year.survR)
    sigma.year.flower ~ dunif(0, 5)
    tau.year.flower <- 1 / (sigma.year.flower * sigma.year.flower)
    sigma.year.germ ~ dunif(0, 5)
    tau.year.germ <- 1 / (sigma.year.germ * sigma.year.germ)

    # density dependenct vital rates with temp/precip effects of year t and t-1 
    # varying slopes by site
    for (v in 1:15){
      beta.mean[v] ~ dunif(-5, 5)
      sigma.beta[v] ~ dunif(0, 5)
      tau.beta[v] <- 1 / (sigma.beta[v] * sigma.beta[v])
      for (i in 1:nSites){
       site.beta[i, v] ~ dnorm(beta.mean[v], tau.beta[v])
      }
    }

    # Vital rates (annual and site variation, no density dependence)
    for (i in 1:nSites){
       # site.survR[i] ~ dnorm(mean.survR, tau.site.survR)
       # site.flower[i] ~ dnorm(mean.flower, tau.site.flower)
       # site.germ[i] ~ dnorm(mean.germ, tau.site.germ)
        
      # flowering every year, other vital rates nYears-1
      logit(flower[i,1]) <- norm.flower[i,1]
       norm.flower[i,1] ~ dnorm(mean.flower, tau.year.flower)

      # meanrosette[i] <- sum(R[i,]) / nYears

      for (t in 2:nYears){  
          # ddrosette[i, t-1] <- R[i,t-1] - meanrosette[i]
         # ddsurvR[i,t-1] <- survR[i,t-1] / (1 + site.beta[i,16] * R[i,t-1])
         logit(survR[i,t-1]) <- norm.survR[i,t-1]
         norm.survR[i,t-1] ~ dnorm(link.survR[i,t-1], tau.year.survR)T(0.85,2.2) # truncated to keep to realistic values
         link.survR[i,t-1] <- mean.survR +        
                                            site.beta[i,1] * zprec[t-1] + 
                                            site.beta[i,2] * zprec[t] + 
                                            site.beta[i,3] * ztemp[t-1] + 
                                            site.beta[i,4] * ztemp[t] +
                                            site.beta[i,5] * R[i,t-1]
        
        # ddflower[i,t] <- flower[i,t] / (1 + site.beta[i,17] * R[i,t-1])
        logit(flower[i,t]) <- norm.flower[i,t] 
        norm.flower[i,t] ~ dnorm(link.flower[i,t], tau.year.flower)T(-5,-1)
        link.flower[i,t]  <- mean.flower +        
                                            site.beta[i,6] * zprec[t-1] + 
                                            site.beta[i,7] * zprec[t] + 
                                            site.beta[i,8] * ztemp[t-1] + 
                                            site.beta[i,9] * ztemp[t] +
                                            site.beta[i,10] * R[i,t-1]


    
        # ddgerm[i,t-1] <- germ[i,t-1] / (1 + site.beta[i,18] * R[i,t-1])
        log(germ[i,t-1]) <- norm.germ[i,t-1]
        norm.germ[i,t-1] ~ dnorm(link.germ[i,t-1], tau.year.germ)T(-3,3.5)
        link.germ[i,t-1] <- mean.germ +           
                                            site.beta[i,11] * zprec[t-1] + 
                                            site.beta[i,12] * zprec[t] + 
                                            site.beta[i,13] * ztemp[t-1] + 
                                            site.beta[i,14] * ztemp[t] +
                                            site.beta[i,15] * R[i,t-1]


        }
      }

    for(i in 1:nSites) {
      # starting population size
      startN[i] ~ dnegbin(P, alpha)
      R[i,1] <- startN[i] - Fl[i,1]
      Fl[i,1] ~ dbin(flower[i,1], startN[i])

      # process model of demographic transitions
      for(t in 2:nYears) {
        Rsurv[i,t-1] ~ dbin(survR[i,t-1], R[i,t-1])
        # Fsurv[i,t-1] ~ dbin(survF[i,t-1], Fl[i,t-1]) # very small
        Fl[i,t] ~ dbin(flower[i,t], Rsurv[i,t-1])
        FtoR[i,t-1] ~ dpois(Fl[i,t-1]*germ[i,t-1] + mean.iota)
  
        R[i, t] <- FtoR[i,t-1] + Rsurv[i,t-1] - Fl[i,t]

      }
      # flowering ratio and observation of flowers
      for(t in 1:nYears) {
        # y[i,t] ~ dbin(0.89, Fl[i,t])
        y[i,t] ~ dpois(Fl[i,t])
        }
    }
    
    # derived predictions summed across creeks
    for (t in 1:nYears){
      creek[,t] <- c(sum(N[c(1:7), t]), sum(N[c(8:12), t]), sum(N[c(13:14), t]))
      total[t] <- sum(creek[,t])
      for (i in 1:nSites){
        N[i, t] <- R[i, t] + Fl[i, t]
        }
      }
    }
    ", fill=TRUE)
sink()



dat <- list(nSites=nSites, nYears=nYears, y=ydata, ztemp = ztemp, zprec = zprec, 
            creek = matrix(data = NA, nrow = 3, ncol = nYears),
            total = rep(NA, nYears))

pars <- c("lambda", "alpha", "N", "survR", "germ", "flower", "mean.iota", 
          "mean.flower", "sigma.year.flower", 
          "mean.survR", "sigma.year.survR", 
          "mean.germ", "sigma.year.germ", 
          "beta.mean", "sigma.beta", "site.beta",
          "creek", "total")

pars <- c("lambda", "alpha", "mean.iota", 
          "mean.flower", "sigma.year.flower", 
          "mean.survR", "sigma.year.survR", 
          "mean.germ", "sigma.year.germ", 
          "beta.mean", "sigma.beta")


jm <- jags.model("gaura.ipm1.txt", dat,  n.chains=3, n.adapt=250000)
update(jm, n.iter=100000)
jc <- coda.samples(jm, pars, n.iter=10000, thin = 10)
saveRDS(jc, "ipm1_mcmc.rds")


# check MCMC ----
summary(jc)
windows(record=T)
plot(jc)
warnings()
gelman.diag(jc)

g <- matrix(NA, nrow=nvar(jc), ncol=2)
for (v in 1:nvar(jc)) {
    g[v,] <- gelman.diag(jc[,v])$psrf
    }


# Dotplot of covariate effects
MCMCplot(jc, params = 'beta.mean',ref_ovl = TRUE,
         labels = c(paste("Rosette survival:", c('Precip (t-1)', 'Precip (t)', 
                    'Temp (t-1)', 'Temp (t)', 'N (t-1)'), sep = " "),
                    paste("Flowering rate:", c('Precip (t-1)', 'Precip (t)', 
                                                'Temp (t-1)', 'Temp (t)', 'N (t-1)'), sep = " "),
                    paste("Germination:", c('Precip (t-1)', 'Precip (t)', 
                                                'Temp (t-1)', 'Temp (t)', 'N (t-1)'), sep = " ")))   #[c(1, 7, 13, 2, 8, 14, 3, 9, 15, 4, 10, 16, 5, 11, 17, 6, 12, 18)])


# simulate future populations and extinction risk ----
nsim <- 500
nyr <- 12

# different conditions in 4 scenarios: random year, random within each of 3 decades
weather2018 <- data.frame(year = year1:year3,
                          temp = ztemp,
                          prec = zprec) %>% 
  filter(year == 2018) 

weather <- data.frame(year = year1:year3,
                      temp = ztemp,
                      prec = zprec) %>% 
  filter(year < 2019)# change years of weather to draw from
# filter(year < 2009 & year > 1998)# drought years

# weather draws
ntemp <- array(sample(x = weather$temp, size = nsim * nyr, replace = TRUE), dim = c(nsim, nyr))
nprec <- array(sample(x = weather$prec, size = nsim * nyr, replace = TRUE), dim = c(nsim, nyr))


# posterior population/parameter estimates
ex <- MCMCchains(jc)
ex <- MCMCchains(jc)[sample(nrow(ex), size = nsim, replace = FALSE),]

# for each draw

# last year population
ns <- ex[, grepl(colnames(ex), pattern = "N[", fixed = TRUE) & grepl(colnames(ex), pattern = ",30]", fixed = TRUE)]
site.r <- ex[, grepl(colnames(ex), pattern = "site.r[", fixed = TRUE)]
site.b <- ex[, grepl(colnames(ex), pattern = "site.beta[", fixed = TRUE)]
iota <- ex[, grepl(colnames(ex), pattern = "iota", fixed = TRUE)]
mean.flower <- ex[, grepl(colnames(ex), pattern = "mean.flower", fixed = TRUE)]
mean.survR <- ex[, grepl(colnames(ex), pattern = "mean.survR", fixed = TRUE)]
mean.germ <- ex[, grepl(colnames(ex), pattern = "mean.germ", fixed = TRUE)]
sigma.year.survR <- ex[, grepl(colnames(ex), pattern = "sigma.year.survR", fixed = TRUE)]
sigma.year.flower <- ex[, grepl(colnames(ex), pattern = "sigma.year.flower", fixed = TRUE)]
sigma.year.germ <- ex[, grepl(colnames(ex), pattern = "sigma.year.germ", fixed = TRUE)]


N <- Fl <- R <- array(NA, dim = c(nrow(ns), ncol(ns), nyr))
for (j in 1:nrow(ns)){  #MCMC draw
  for (i in 1:ncol(ns)){  # site
    N[j,i,1] <- ns[j,i]
    ntemp[j,1] <- weather2018$temp
    nprec[j,1] <- weather2018$prec
    
    Fl[j,i,1] <-  rbinom(1, N[j,i,1], plogis(rnorm(1, mean.flower[j], sigma.year.flower[j])))
    R[j,i,1] <- N[j,i,1] - Fl[j,i,1]
    
    
    for (t in 2:nyr){ # year
      
      link.survR <- mean.survR[j] +  
        site.b[j, c(1:14 + 14 * (1-1))[i]] * nprec[j,t-1] + 
        site.b[j, c(1:14 + 14 * (2-1))[i]] * nprec[j,t] + 
        site.b[j, c(1:14 + 14 * (3-1))[i]] * ntemp[j,t-1] + 
        site.b[j, c(1:14 + 14 * (4-1))[i]] * ntemp[j,t] +
        site.b[j, c(1:14 + 14 * (5-1))[i]] * R[j,i,t-1]
      norm.survR <- rnorm(1, link.survR, sigma.year.survR[j])
      survR <- plogis(norm.survR)
      
      link.flower  <- mean.flower[j] +  
        site.b[j, c(1:14 + 14 * (6-1))[i]] * nprec[j,t-1] + 
        site.b[j, c(1:14 + 14 * (7-1))[i]] * nprec[j,t] + 
        site.b[j, c(1:14 + 14 * (8-1))[i]] * ntemp[j,t-1] + 
        site.b[j, c(1:14 + 14 * (9-1))[i]] * ntemp[j,t] +
        site.b[j, c(1:14 + 14 * (10-1))[i]] * R[j,i,t-1]
      norm.flower <- rnorm(1, link.flower, sigma.year.flower[j])
      flower <- plogis(norm.flower) 
      
      
      
      
      link.germ <- mean.germ[j] +  
        site.b[j, c(1:14 + 14 * (11-1))[i]] * nprec[j,t-1] + 
        site.b[j, c(1:14 + 14 * (12-1))[i]] * nprec[j,t] + 
        site.b[j, c(1:14 + 14 * (13-1))[i]] * ntemp[j,t-1] + 
        site.b[j, c(1:14 + 14 * (14-1))[i]] * ntemp[j,t] +
        site.b[j, c(1:14 + 14 * (15-1))[i]] * R[j,i,t-1]
      norm.germ <- rnorm(1, link.germ, sigma.year.germ[j])
      germ <- exp(norm.germ)
      
      Rsurv <- rbinom(1, R[j,i,t-1], survR)
      Fl[j,i,t] <- rbinom(1, Rsurv, flower)
      FtoR <- rpois(1, Fl[j,i,t-1]*germ + iota[j])
      
      R[j,i, t] <- FtoR + Rsurv - Fl[j,i,t]
      N[j,i,t] <- R[j,i,t] + Fl[j,i,t]
      
    }
  }
}
# N <- Fl
# by segment
forecast <- apply(N, MARGIN = c(2,3), quantile, probs = c(0.05, 0.1, 0.25, .5, 0.75, 0.9, 0.95), na.rm = TRUE) 
extinct <- apply(N, MARGIN = c(2,3), function(x) length(which(x == 0))/length(x)) 

segs <- census %>%
  filter(segment %!in% c("C-VIII")) %>% 
  droplevels()
segnames <- levels(segs$segment)


flist <- list()
elist <- list()
for(y in 1:dim(forecast)[3]){
  ftmp <- data.frame(t(forecast[,,y])) %>% 
    mutate_all(as.numeric) %>% 
    mutate(segment = segnames,
           year = 2018 + y)
  flist[[y]] <- ftmp
  
  etmp <- data.frame(segment = segnames, year = 2018 + y, prob_ext = extinct[,y])
  elist[[y]] <- etmp
}
fdf <- bind_rows(flist)
edf <- bind_rows(elist)


# by creek
CC <- abind::abind(apply(N[,1:7,], MARGIN = c(1,3), sum))
DC <- abind::abind(apply(N[,8:12,], MARGIN = c(1,3), sum))
UC <- abind::abind(apply(N[,13:14,], MARGIN = c(1,3), sum))
# total
ALL <- abind::abind(apply(N[,1:14,], MARGIN = c(1,3), sum))
res <- list(CC, DC, UC, ALL)
resname <- c("Crow", "Diamond", "Unnamed", "All")
flist <- list()
elist <- list()
for (i in 1:4){
  
  ftmp <- data.frame(t(apply(res[[i]], MARGIN = 2, quantile, probs = c(0.05, 0.1, 0.25, .5, 0.75, 0.9, 0.95), na.rm = TRUE))) %>% 
    mutate_all(as.numeric) %>% 
    mutate(segment = resname[i],
           year = 2019:2030)
  flist[[i]] <- ftmp
  extinct <- data.frame(year = 2019:2030, segment = resname[i],
                        prob_ext = apply(res[[i]], MARGIN = 2, function(x) length(which(x == 0))/length(x)))
  elist[[i]] <- extinct
}
summ_fdf <- bind_rows(fdf, bind_rows(flist))
summ_edf <- bind_rows(edf, bind_rows(elist))

summ_fdf <- bind_rows(flist)

# just plot summarized by creek
plt <- ggplot(data = summ_fdf, aes(x = year, y = X50.)) +
  geom_point(size = 2) +
  # geom_linerange(aes(ymin = X5., ymax = X95.), size = 1, alpha = .25) +
  # geom_linerange(aes(ymin = X25., ymax = X75.), size = 1.5, alpha = .5) +
  facet_wrap(~segment, scales = "free", ncol = 2) +
  ylab("Forecasted population size") +
  ggtitle("Forecasted population size with 1989-2018 weather")
plt  


out_edf <- summ_edf %>% 
  filter(year == 2030)
write.csv(out_edf, row.names = FALSE, file = "ext_flower_normal.csv")








# plot uncertainty ----

# could also grab rows from mcmclist, sort by variables, and plot population trajectories over time

ch <- rbind(jc[[1]], jc[[2]], jc[[3]])
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

ks <- df[grep(df$vars, pattern = "K[", fixed = TRUE),] %>% 
  mutate(segment = segnames)

creek <- df[grep(df$vars, pattern = "creek", fixed = TRUE),] %>% 
  mutate(segment = rep(c("crow", "diamond", "unnamed"), nYears),
         year = 1989 - 1 + as.numeric(gsub("\\D+", "", str_split_fixed(vars, ",", 2)[,2])))

total <- df[grep(df$vars, pattern = "total", fixed = TRUE),] %>% 
  mutate(segment = c("all"),
         year = 1989:year3 )

ns <- bind_rows(ns, creek, total)


# ps <- df[grep(df$vars, pattern = "p.mean", fixed = TRUE),] %>% 
#   mutate(segment = as.numeric(gsub("\\D+", "", str_split_fixed(vars, ",", 2)[,1])),
#          year = 1989 - 1 + as.numeric(gsub("\\D+", "", str_split_fixed(vars, ",", 2)[,2])))
# 
# ext <- df[grep(df$vars, pattern = "PQE[", fixed = TRUE),] %>% 
#   mutate(segment = as.numeric(gsub("\\D+", "", str_split_fixed(vars, ",", 2)[,1])),
#          year = 1989 - 1 + as.numeric(gsub("\\D+", "", str_split_fixed(vars, ",", 2)[,2]))) %>% 
#   filter(segment == 3)

# nus <- df[grep(df$vars, pattern = "nu[", fixed = TRUE),]
fs <- df[grep(df$vars, pattern = "flower[", fixed = TRUE),]

rs <- df[grep(df$vars, pattern = "R[", fixed = TRUE),][1:448,]

bs <- df[grep(df$vars, pattern = "beta", fixed = TRUE),]

sr <- df[grep(df$vars, pattern = "survR[", fixed = TRUE),]

sg <- df[grep(df$vars, pattern = "sigma", fixed = TRUE),]

ms <- df[grep(df$vars, pattern = "mean", fixed = TRUE),]


plt <- ggplot(data = ns %>% filter(year < 2026), aes(x = year, y = X50.)) +
  geom_point() +
  geom_linerange(aes(ymin = X25., ymax = X75.), alpha = 0.5) +
  geom_point(data = counts %>% filter(year < 2026), aes(x = year, y = flower.count), color = "red") +
  # geom_hline(data = ks, aes(yintercept = X50.), linetype = "dashed") +
  # geom_hline(data = ks, aes(yintercept = X25.), linetype = "dashed", color = "grey") +
  # geom_hline(data = ks, aes(yintercept = X75.), linetype = "dashed", color = "grey") +
  facet_wrap(~segment, scales = "free", ncol = 3) +
  ylab("Count") +
  ggtitle("Median prediction of total population size (with 50% CI) versus flower census")
plt  


test <- left_join(ns, ydata, by = c("year", "segment"))

plot(nus$X50., ps$X50.)

