library(lme4)
packageVersion("lme4")  ## 1.1.4, but this should work as long as >1.0.0
library(MuMIn)

cbpp <- transform(cbpp,prop=incidence/size)

gm0 <- glmer(cbind(incidence, size - incidence) ~ period+(1|herd),
             family = binomial, data = cbpp)
dredge.models<-dredge(gm0,trace=FALSE,rank="AICc")
my.dredge.models<-get.models(dredge.models)
silly<-model.avg(my.dredge.models,subset=delta<10)

predict(silly,type="response")

library(ggplot2)
theme_set(theme_bw())  ## cosmetic
g0 <- ggplot(cbpp,aes(period,prop))+
  geom_point(alpha=0.5,aes(size=size))

predframe <- data.frame(period=levels(cbpp$period))

predframe$prop <- predict(gm0,newdata=predframe,type="response",ReForm=NA)

g0 + geom_point(data=predframe,colour="red")+
  geom_line(data=predframe,colour="red",aes(group=1))







setwd("C:/Users/Tyson/Dropbox/Gaura")
data <- read.csv("veg.csv", header = TRUE)

data$id <- paste(data$Creek, data$Subunit, sep = ".")
data$id <- as.factor(data$id)

data$prop.flower <- data$Flow04 / (data$Veg04 + data$Flow04)
data$total <- data$Flow04 + data$Veg04

#1st attempt, not right
#mod <- glm(Flow04 ~ total + id, family = poisson, data = data)

#code adapted from GLMM from Ap Browns
surv.var =cbind(data$Flow04, data$total - data$Flow04)

#random effects probably not necessary, has lower AIC, but extra complication not changing parameters
mod <- glmer(surv.var ~ 1 + (1|id), data = data, 
             family = binomial, control=glmerControl(optimizer="bobyqa"))

mod1 <- glm(surv.var ~ 1 + Creek + log(total), data = data, 
            family = binomial)

predict(mod,type="response")

library(ggplot2)
theme_set(theme_bw())  ## cosmetic
g0 <- ggplot(data,aes(Creek, prop.flower))+
  geom_point(alpha=0.5,aes(size=log(total))) 
g0 + labs(x = "Creek", y = "Percent of plants flowering in plot")


tot.pred <- seq(1,200, 1)

predframeCrow <- data.frame(total = tot.pred, Creek = factor(rep("Crow", length(tot.pred)), levels = levels(data$Creek)))
predframeCrow$prop.flower <- predict(mod1,newdata=predframeCrow,type="response",ReForm=NA)

predframeDiam <- data.frame(total = tot.pred, Creek = factor(rep("Diamond", length(tot.pred)), levels = levels(data$Creek)))
predframeDiam$prop.flower <- predict(mod1,newdata=predframeDiam,type="response",ReForm=NA)

predframeUnna <- data.frame(total = tot.pred, Creek = factor(rep("Unnamed", length(tot.pred)), levels = levels(data$Creek)))
predframeUnna$prop.flower <- predict(mod1,newdata=predframeUnna,type="response",ReForm=NA)


plot(c(1,200), c(0,1), type = "n", xlab = "#plants",
     ylab = "prop flowering", log = "x")
lines(predframeCrow$prop.flower, col = "red")
lines(predframeDiam$prop.flower, col= "blue")
lines(predframeUnna$prop.flower, col = "green")


library(ggplot2)
hist_cut <- ggplot(data, aes(x=prop.flower, fill = Creek))
hist_cut + geom_bar(position="dodge")

library(boot)
#quantify variation in proportions
(coeftbl <- as.data.frame(coef(summary(mod))))
## 95% confidence intervals
out <- with(coeftbl, Estimate + outer(`Std. Error`, c(lower=-1, upper=1)) * sqrt(qchisq(0.95, 1)))
inv.logit(mod@beta)
inv.logit(out)


#Does the Floyd and Ranker data match the pattern?







