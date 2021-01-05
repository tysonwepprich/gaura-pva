#gaura matrix model using floyd & ranker raw data
library(rjags)
library(R2jags)
library(popbio)
library(plyr)
library(reshape)
library(lme4)

setwd("C:/Users/Tyson/Dropbox/Gaura")
data <- read.csv("matrix.data", header = TRUE)
flowers <- read.csv("floyd.flowers.csv", header = TRUE)
flowers$Creek <- c(rep("Crow", 9), rep("Unnamed",9), rep("Diamond", 9))


###################################
#Seed survial and germination rates from new seedlings in demographic plots

seed.data <- ddply(data, .(plot, year.start), summarize,
                   seedlings = sum(num.plants[stage == "seed"]),
                   flowers = sum(num.plants[stage == "flower"]),
                   total.plants = sum(num.plants))
seed.data$plot <- as.factor(seed.data$plot)

seed.data$Creek <- c(rep("Crow", 6), rep("Unnamed",6), rep("Diamond", 6))

mod <- lm(seedlings ~ flowers + total.plants + as.factor(year.start), data = seed.data)
mod2 <- lm(seedlings[year.start == 1993] ~ flowers[year.start == 1992] + flowers[year.start == 1993], data = seed.data)
mod.mix <- lmer(seedlings ~ flowers + total.plants + as.factor(year.start) + (1|plot), data = seed.data)

surv.var =cbind(seed.data$flowers, seed.data$total.plants - seed.data$flowers)

mod <- glmer(surv.var ~ 1 + Creek + log(total.plants) + factor(year.start) + (1|plot), data = seed.data, 
             family = binomial, control=glmerControl(optimizer="bobyqa"))

mod1 <- glm(surv.var ~ 1 + Creek + log(total.plants) + factor(year.start), data = seed.data, 
            family = binomial)

dredge.models<-dredge(mod.mix,trace=FALSE,rank="AICc")
my.dredge.models<-get.models(dredge.models)
silly<-model.avg(my.dredge.models,subset=delta<10)

flowers$prop <- flowers$flowers/flowers$total
#using flowers data
prop.flow <- cbind(flowers$flowers, flowers$total - flowers$flowers)
mod <- glm(prop.flow ~ 1 + log(total) + factor(year) + Creek, data = flowers, family = binomial)
mod.mix <- glmer(prop.flow ~ 1 + log(total) + factor(year) + (1|plot), data=flowers,
                 family = binomial, control=glmerControl(optimizer="bobyqa"))

#quantify variation in proportions
mod.var <- glmer(prop.flow ~ 1 + (1|plot), data = flowers, family = binomial, control=glmerControl(optimizer="bobyqa"))
(coeftbl <- as.data.frame(coef(summary(mod.var))))
## 95% confidence intervals
out <- with(coeftbl, Estimate + outer(`Std. Error`, c(lower=-1, upper=1)) * sqrt(qchisq(0.95, 1)))
inv.logit(mod.var@beta)
inv.logit(out)


library(ggplot2)
hist_cut <- ggplot(flowers, aes(x=prop, fill = Creek))
hist_cut + geom_bar(position="dodge")

flowers$year <- as.factor(flowers$year)
theme_set(theme_bw())  ## cosmetic
g0 <- ggplot(flowers,aes(year, prop))+
  geom_point(alpha=0.5,aes(size=log(total)))
g0 + labs(x = "Year", y = "Percent of plants flowering in plot")
library(boot)
inv.logit(mod.mix@beta)
###################################


expanded<-untable(data[, c(2:5)], data[,1])
expanded$stage <- factor(expanded$stage, levels = c("seed", "small", "medium", "large1", "large2", "flower", "dead"), ordered = TRUE)
expanded$fate <- factor(expanded$fate, levels = c("seed", "small", "medium", "large1", "large2", "flower", "dead"), ordered = TRUE)

stages <- c("seed", "small", "medium", "large1", "large2", "flower")
fates <- c("dead", "small", "medium", "large1", "large2", "flower")

tf <- table(expanded[, "fate"], expanded[, "stage"])
T.mat <- prop.table(tf, 2)[fates, stages]


all.mats <- array(data = NA, dim = c(9, 2, 6, 6))
for (plot in 1:9){
  for (year in 1992:1993){
    temp <- expanded[expanded$plot == plot & expanded$year.start == year, ]
    tf <- table(temp[, "fate"], temp[, "stage"])
    year.index <- year - 1991
    all.mats[plot, year.index, ,] <- prop.table(tf, 2)[fates, stages]
  }
}


#using popbio functions
test <- expanded
test <- expanded[expanded$plot == 5 & expanded$year.start == 1992, ]
test$fruits <- 0
test$fruits[test$stage == "flower"] <- rpois(length(test$stage[test$stage == "flower"]), 50)
recruits  <- length(test$stage[test$stage == "seed"])

seed.bank.size = 1000
seed.survival = .5
seeds.per.fruit = 3
seeds.from.plants <- sum(test$fruits) * seeds.per.fruit
recruitment.rate <- recruits/(seed.bank.size + seeds.from.plants) 
test$recruit <- test$fruits/sum(test$fruits) * seeds.from.plants * recruitment.rate 
test$seed <- test$fruits * seeds.per.fruit * seed.survival

mod <- projection.matrix(test)
mod[1,1] <- seed.survival * (1 - recruitment.rate)
mod[2:5,1] <- mod[2:5,1] * recruitment.rate
mod[7,1] <- 1 - seed.survival
mod[7,7] <- 1

A <- mod
lam <- max(Re(eigen(A)$values))
n <- summary(test$stage)
n1 <- A %*% n
n1 <- as.vector(n1)
n2 <- as.vector(A %*% n1)
n3 <- as.vector(A %*% n2)



#TRY it out in JAGS

#specify model in JAGS language
sink("matrix.model.txt")
cat("
    model{
    
    #Parameters
    #surv: Survival for each stage
    #psi: transitions from each stage
    #fert: fertility/flower
    #rec: recruitment rate from seed
    #sb: size of seedbank

    #States:
#     1 seed
#     2 small
#     3 medium
#     4 large1
#     5 large2
#     6 flower
#     7 dead
    
    #Priors and constraints
      #Survival for each stage
      for (i in 1:6){
        surv[i] ~ dunif(0,1)
      }
    
      #Transitions
      for (i in 1:5){
        a[i] ~ dgamma(1,1)
        psi.sd  <- a[i]/sum(a[])
        b[i] ~ dgamma(1,1)
        psi.sm  <- b[i]/sum(b[])
        c[i] ~ dgamma(1,1)
        psi.m  <- c[i]/sum(c[])
        d[i] ~ dgamma(1,1)
        psi.l1  <- d[i]/sum(d[])
        e[i] ~ dgamma(1,1)
        psi.l2  <- e[i]/sum(e[])
        f[i] ~ dgamma(1,1)
        psi.fl <- f[i]/sum(f[])
      }

      #Other parameters
      fert ~ dunif(0, 1000)
      rec ~ dunif(0,1)
      sb ~ dnorm(100, 0.001)T(0,)

    #Define probability of state S(t+1) given S(t) for multinomial draws from observed age-structured plants
      mat[1,1] <- surv[1] * (1 - rec)
      mat[1,2] <- 0
      mat[1,3] <- 0
      mat[1,4] <- 0
      mat[1,5] <- 0
      mat[1,6] <- fert * surv[1]
      for (i in 2:6){
        mat[i,1] <- rec * surv[1] * psi.sd[i-1]
      }
      for (i in 2:6){
        mat[i,2] <- surv[2] * psi.sm[i-1]
      }
      for (i in 2:6){
        mat[i,3] <- surv[3] * psi.m[i-1]
      }
      for (i in 2:6){
        mat[i,4] <- surv[4] * psi.l1[i-1]
      }      
      for (i in 2:6){
        mat[i,5] <- surv[5] * psi.l2[i-1]
      }      
      for (i in 2:6){
        mat[i,6] <- surv[6] * psi.fl[i-1]
      }      
      for (i in 1:6){
        mat[7,i] <- 1 - surv[i]
        mat[i,7] <- 0
      }
      mat[2,5] <- 0
      mat[2,6] <- 0
      mat[3,6] <- 0
      mat[6,1] <- 0
      mat[7,7] <- 1

    #Likelihood
    



    }")

