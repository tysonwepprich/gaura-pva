# vegetative plant frequency

setwd("C:/Users/Tyson/Dropbox/Gaura")
data <- read.csv("veg.csv", header = TRUE)

data$id <- paste(data$Creek, data$Subunit, sep = ".")
data$id <- as.factor(data$id)

data$prop.flower <- data$Flow04 / (data$Veg04 + data$Flow04)
data$total <- data$Flow04 + data$Veg04

mod <- glm(Flow04 ~ total + id, family = poisson, data = data)

#code adapted from GLMM from Ap Browns
surv.var =cbind(data$Flow04, data$total - data$Flow04)


mod <- glmer(surv.var ~ 1 + log(total) + (1|id), data = data, 
             family = binomial, control=glmerControl(optimizer="bobyqa"))

mod <- glmer(surv.var ~ 1 + (1|Creek), data = data, family = binomial)


library(ggplot2)
hist_cut <- ggplot(data, aes(x=prop.flower, fill = Creek))
hist_cut + geom_bar(position="dodge")

library(boot)
inv.logit(mod.mix@beta)




gaura <- read.csv("subcounts.csv", header = TRUE)
gaura.id <- as.factor(paste(gaura$Creek, gaura$Subunit, sep = "."))
pop.names <- levels(reorder(levels(gaura.id), c(8, 9, 10, 11, 12, 13, 14, 3, 4, 5, 6, 7, 1, 2)))

gaura.mat <- matrix(data = gaura$Count, nrow = 25)
colnames(gaura.mat) <- pop.names
rownames(gaura.mat) <- c(1989:2013)

gaura.dat <- log(t(gaura.mat + 1))

#MARSS

library(MARSS)

Z.model <- factor(c(rep("U",2), rep("D",5), rep("C",7)))
Q.model <- "unconstrained"
B.model <- "diagonal and unequal"

fit <- MARSS(gaura.dat, model = list(Q = Q.model, Z = Z.model, B = B.model))

#First, test if structure between 3 subpopulations
unnamed <- log(rowSums(gaura.mat[, c(1:2)], na.rm = TRUE))
diamond <- log(rowSums(gaura.mat[, c(3:7)], na.rm = TRUE))
crow <- log(rowSums(gaura.mat[, c(7:14)], na.rm = TRUE))

plot(1989:2013, crow, type = "l", col = "red", xlim = c(1989, 2013), ylim = c(2, 9))
lines(1989:2013, diamond, col = "blue")
lines(1989:2013, unnamed, col = "green")

pops <- cbind(unnamed, diamond, crow)

dat <- as.matrix(t(pops))

Z1 <- factor(c(1,1,1))
Z2 <- factor(c(1:3))

Z.models <- list(Z1, Z2)
names(Z.models) <- c("one pop", "three subpop")
#three subpopulations makes more sense
Q.models <- c("diagonal and equal","diagonal and unequal", "unconstrained")
#unconstrained to see correlations of process error
U.models <- c("equal","unequal", "zero")
#makes sense for trend to be unequal
R.models <- c("diagonal and equal", "diagonal and unequal")
#makes sense for R (obs error) to be equal
A.model <- "scaling"
B.models <- c("identity", "diagonal and equal", "diagonal and unequal")
#makes sense to be diagonal and unequal (density dependence)
x0.model <- "unequal"
V0.model <- "zero"
model.constant=list(A=A.model, x0=x0.model, V0=V0.model, tinitx=0)

#data mining to try it out
out.tab=NULL
fits=list()
for(i in 1:length(Z.models)){
  for(Q.model in Q.models){
    for(U.model in U.models){
      for(R.model in R.models){
        for(B.model in B.models){
    fit.model = c(list(Z=Z.models[[i]], Q=Q.model, U=U.model, R=R.model, B=B.model), model.constant)
    fit = MARSS(dat, model=fit.model,
                silent=TRUE, control=list(maxit=1000))
    out=data.frame(H=names(Z.models)[i], Q=Q.model, U=U.model, R=R.model, B=B.model,
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   m=length(unique(Z.models[[i]])),
                   num.iter=fit$numIter, converged=!fit$convergence)
    out.tab=rbind(out.tab,out)
    fits=c(fits,list(fit))
    if(i==5) next #one m for panmictic so only run 1 Q
  }
}
}}}

min.AICc=order(out.tab$AICc)
out.tab.1=out.tab[min.AICc,]
out.tab.1=cbind(out.tab.1, delta.AICc=out.tab.1$AICc-out.tab.1$AICc[1])
out.tab.1=cbind(out.tab.1, rel.like=exp(-1*out.tab.1$delta.AICc/2))
out.tab.1=cbind(out.tab.1, AIC.weight = out.tab.1$rel.like/sum(out.tab.1$rel.like))



#do it again with more constraints on parameters

Z.model <- factor(c(1:3))
#three subpopulations makes more sense
Q.model <- "diagonal and equal"
#unconstrained to see correlations of process error
U.models <- c("unequal", "equal", "zero")
#makes sense for trend to be unequal
R.model <- "diagonal and equal"
#makes sense for R (obs error) to be equal
A.model <- "equal"
#having lower A for herbivory years better
B.model <- c("identity")

#C.models <- c("unconstrained", "zero")
#covs <- list(matrix(covariates[1,], nrow = 1), matrix(covariates[2,], nrow=1), matrix(covariates[3,], nrow=1), matrix(covariates[4,], nrow=1),
#             covariates[1:2,], covariates[2:3,], covariates[3:4,])

#makes sense to be diagonal and unequal (density dependence)
x0.model <- "unequal"
V0.model <- "zero"
model.constant=list(B=B.model, Z=Z.model, R=R.model, A=A.model, Q=Q.model, x0=x0.model, V0=V0.model, tinitx=1)

#data mining to try it out
out.tab=NULL
fits=list()
        for(U.model in U.models){
          fit.model = c(list(U=U.model), model.constant)
          fit = MARSS(dat, model=fit.model,
                      silent=TRUE, control=list(maxit=2000, safe = TRUE, trace = 1))
          out=data.frame(B=B.model,
                         logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                         num.iter=fit$numIter, converged=!fit$convergence)
          out.tab=rbind(out.tab,out)
          fits=c(fits,list(fit))
        }
      
    


min.AICc=order(out.tab$AICc)
out.tab.1=out.tab[min.AICc,]
out.tab.1=cbind(out.tab.1, delta.AICc=out.tab.1$AICc-out.tab.1$AICc[1])
out.tab.1=cbind(out.tab.1, rel.like=exp(-1*out.tab.1$delta.AICc/2))
out.tab.1=cbind(out.tab.1, AIC.weight = out.tab.1$rel.like/sum(out.tab.1$rel.like))



plot(1989:2013, fits[[1]]$states[1,], type = "l", col = "red", xlim = c(1989, 2013), ylim = c(2, 9))
lines(1989:2013, fits[[1]]$states[2,], col = "blue")
lines(1989:2013, fits[[1]]$states[3,], col = "green")

CIs.1=MARSSparamCIs(fits[[1]])
CIs.2=MARSSparamCIs(fits[[2]])

plot.ts(residuals(fits[[2]])$model.residuals[1,],
        ylab="Residual")
abline(h=0, lty="dashed")

plot(covariates[2,], residuals(fits[[2]])$model.residuals[2,])


A.time <- matrix(c("r","r","r"),3,1)
herb <- matrix(c("h", "h", "h"),3,1)
At <- array(A.time, dim = c(dim(A.time), dim(dat)[2]))
At[,,19:20] <- herb

U1=matrix("t1",5,1); U2=matrix("t2",5,1)
Ut=array(U2,dim=c(dim(U1),dim(dat)[2]))
Ut[,,1:floor(TT/2)]=U1





