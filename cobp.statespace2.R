# Gaura State-space
# Setup----
library(dplyr)
library(tidyr)
library(ggplot2)
library(rjags)
# library(runjags)
library(stringr)
# library(dclone)
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



# State-space flowers ----
sink(file="gaura.flower.txt")
cat("
    # Model flower counts as variable of interest (ignoring rosette population)

    model {
    # starting population sizes with negative binomial
    lambda ~ dunif(0, 5000)
    alpha ~ dunif(0, 20)
    P <- alpha/(alpha + lambda)

    iota ~ dunif(0, 10) # immigration
    
    # varying slopes by site for covariates
    for (v in 1:5){
      beta.mean[v] ~ dunif(-5, 5)
      sigma.beta[v] ~ dunif(0, 5)
      tau.beta[v] <- 1 / (sigma.beta[v] * sigma.beta[v])
      for (i in 1:nSites){
        site.beta[i, v] ~ dnorm(beta.mean[v], tau.beta[v])
        }
      }
    
    # growth rate varying by site
    mean.r ~ dnorm(0, 0.001)
    sigma.r ~ dunif(0, 5)
    tau.r <- 1 / (sigma.r * sigma.r)
    for(i in 1:nSites){
      site.r[i] ~ dnorm(mean.r, tau.r)
    }
    
    # environmental stochasticity in growth rates
    sigma.nu ~ dunif(0, 5)
    tau.nu <- 1 / (sigma.nu * sigma.nu)
    for(t in 2:nYears) {
      for(i in 1:nSites){
        nu[i, t-1] ~ dnorm(0, tau.nu)
      }
    }
    
    for(i in 1:nSites) {
      N[i,1] ~ dnegbin(P, alpha)
      for(t in 2:nYears) {
        N[i,t] ~ dpois(muN[i, t-1])
        muN[i,t-1] <- N[i, t-1] * exp(r[i,t-1] + nu[i,t-1]) + iota
    
        r[i,t-1] <-        site.r[i] +     
                                        site.beta[i,1] * zprec[t-1] + 
                                        site.beta[i,2] * zprec[t] + 
                                        site.beta[i,3] * ztemp[t-1] + 
                                        site.beta[i,4] * ztemp[t] +
                                        site.beta[i,5] * N[i,t-1]
      }
      for(t in 1:nYears) {
        # y[i,t] ~ dbin(0.87, N[i,t])
        y[i,t] ~ dpois(N[i,t])
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

dat.fl <- list(nSites=nSites, nYears=nYears, y=ydata, ztemp = ztemp, zprec = zprec, zyear = zyear,
                   creek = matrix(data = NA, nrow = 3, ncol = nYears),
                   total = rep(NA, nYears))

pars.fl <- c("lambda", "alpha", "r", "sigma.nu", "iota", "N", "site.r",
                 "sigma.beta", "beta.mean", "site.beta", "mean.r", "sigma.r", "creek", "total")

jm.fl <- jags.model("gaura.flower.txt", dat.fl, n.chains=3, n.adapt=50000)

# update(jm.fl, n.iter=50000)
jc.fl <- coda.samples(jm.fl, pars.fl, n.iter=1000, thin = 1)


saveRDS(jc.fl, "flower_mcmc.rds")




# check MCMC ----
summary(jc.fl)
windows(record=T)
plot(jc.fl)
warnings()
gelman.diag(jc.fl)


jc <- jc.fl
# gelman diagnostics for each variable separately
g <- matrix(NA, nrow=nvar(jc), ncol=2)
for (v in 1:nvar(jc)) {
  g[v,] <- gelman.diag(jc[,v])$psrf
}

MCMCsummary(jc, round = 2)

# Dotplot of covariate effects
MCMCplot(jc, params = 'beta.mean',ref_ovl = TRUE,
         labels = c('Precipitation (t-1)', 'Precipitation (t)', 
                    'Temperature (t-1)', 'Temperature(t)', 'Population Size (t-1)'))

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


ex <- MCMCchains(jc)
# posterior population/parameter estimates
ex <- MCMCchains(jc)[sample(nrow(ex), size = nsim, replace = FALSE),]

# for each draw

# last year population
ns <- ex[, grepl(colnames(ex), pattern = "N[", fixed = TRUE) & grepl(colnames(ex), pattern = ",30]", fixed = TRUE)]
site.r <- ex[, grepl(colnames(ex), pattern = "site.r[", fixed = TRUE)]
site.b <- ex[, grepl(colnames(ex), pattern = "site.beta[", fixed = TRUE)]
iota <- ex[, grepl(colnames(ex), pattern = "iota", fixed = TRUE)]
sigma.nu <- ex[, grepl(colnames(ex), pattern = "sigma.nu", fixed = TRUE)]

# weather draws
ntemp <- array(sample(x = weather$temp, size = nsim * nyr, replace = TRUE), dim = c(nrow(ns), nyr))
nprec <- array(sample(x = weather$prec, size = nsim * nyr, replace = TRUE), dim = c(nrow(ns), nyr))


N <- array(NA, dim = c(nrow(ns), ncol(ns), nyr))
for (j in 1:nrow(ns)){  #MCMC draw
  for (i in 1:ncol(ns)){  # site
    N[j,i,1] <- ns[j,i]
    ntemp[j,1] <- weather2018$temp
    nprec[j,1] <- weather2018$prec
    
    for (t in 2:nyr){ # year
      
      nu <-  rnorm(1, 0, sigma.nu[j])
      
      r <-        site.r[j,i] +  
        site.b[j, c(1:14)[i]] * nprec[j,t-1] + 
        site.b[j, c(15:28)[i]] * nprec[j,t] + 
        site.b[j, c(29:42)[i]] * ntemp[j,t-1] + 
        site.b[j, c(43:56)[i]] * ntemp[j,t] +
        site.b[j, c(57:70)[i]] * N[j,i,t-1]
      
      muN <- N[j,i,t-1] * exp(r + nu)
      N[j,i,t] <- rpois(1, muN)
    }
  }
}

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
  geom_linerange(aes(ymin = X5., ymax = X95.), size = 1, alpha = .25) +
  geom_linerange(aes(ymin = X25., ymax = X75.), size = 1.5, alpha = .5) +
  facet_wrap(~segment, scales = "free", ncol = 2) +
  ylab("Forecasted population size") +
  ggtitle("Forecasted population size with 1999-2008 weather")
plt  


out_edf <- summ_edf %>% 
  filter(year == 2030)
write.csv(out_edf, row.names = FALSE, file = "ext_flower_normal.csv")






# plot population estimates ----

# could also grab rows from mcmclist, sort by variables, and plot population trajectories over time

ch <- rbind(jc.fl[[1]], jc.fl[[2]], jc.fl[[3]])
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


# # carrying capacity two ways
# ks <- df[grep(df$vars, pattern = "K[", fixed = TRUE),] %>% 
#   mutate(segment = segnames)
# 
# 
# sk <- df[grep(df$vars, pattern = "site.beta[", fixed = TRUE),]
# sr <- df[grep(df$vars, pattern = "site.r[", fixed = TRUE),] %>% 
#   mutate(beta = sk[71:84, "X50."],
#          segment = segnames,
#          K = -X50. / beta)



creek <- df[grep(df$vars, pattern = "creek", fixed = TRUE),] %>% 
  mutate(segment = rep(c("crow", "diamond", "unnamed"), nYears),
         year = 1989 - 1 + as.numeric(gsub("\\D+", "", str_split_fixed(vars, ",", 2)[,2])))

total <- df[grep(df$vars, pattern = "total", fixed = TRUE),] %>% 
  mutate(segment = c("all"),
         year = 1989:year3 )

ns <- bind_rows(ns, creek, total)

rs <- df[grep(df$vars, pattern = "r[", fixed = TRUE),]

bs <- df[grep(df$vars, pattern = "beta", fixed = TRUE),]



plt <- ggplot(data = ns %>% filter(year < 2026), aes(x = year, y = X50.)) +
  geom_point() +
  geom_linerange(aes(ymin = X25., ymax = X75.), alpha = 0.5) +
  geom_point(data = counts %>% filter(year < 2026), aes(x = year, y = flower.count), color = "red") +
  facet_wrap(~segment, scales = "free", ncol = 3) +
  ylab("Count") +
  ggtitle("Median prediction of total population size (with 50% CI) versus flower census")
plt  

