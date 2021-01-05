#final try at Gaura analysis for update



gaura <- read.csv("subcounts.csv", header = TRUE)
gaura.id <- as.factor(paste(gaura$Creek, gaura$Subunit, sep = "."))
pop.names <- levels(reorder(levels(gaura.id), c(8, 9, 10, 11, 12, 13, 14, 3, 4, 5, 6, 7, 1, 2)))

gaura.mat <- matrix(data = gaura$Count, nrow = 25)
colnames(gaura.mat) <- pop.names
rownames(gaura.mat) <- c(1989:2013)

gaura.dat <- log(t(gaura.mat + 1))


#First, run most basic PVA using default settings in MARSS
#Try 3 subpopulation counts independently

unnamed <- cbind(1989:2013, log(rowSums(gaura.mat[, c(1:2)], na.rm = TRUE)))
diamond <- cbind(1989:2013, log(rowSums(gaura.mat[, c(3:7)], na.rm = TRUE)))
crow <- cbind(1989:2013, log(rowSums(gaura.mat[, c(7:14)], na.rm = TRUE)))

#run three times, saved figueres and model fit for report
out <- CSEGriskfigure(crow, te = 50, absolutethresh = FALSE, threshold = .1, datalogged = TRUE, 
                      return.model = TRUE, CI.method = "hessian")


out <- CSEGmod(diamond, te = 50, absolutethresh = FALSE, threshold = .1, datalogged = TRUE, 
               return.model = TRUE, CI.method = "hessian")



#Incorporate density dependence for 3 subpopulations run individually.
dat <- crow[,2]
kem = MARSS(dat, model = list(U = matrix("U"), B = "unconstrained"), silent = TRUE, control=list(maxit=2000))
MARSSparamCIs(kem, method = "hessian")

dat <- diamond[,2]
kem = MARSS(dat, model = list(U = matrix("U"), B = "unconstrained"), silent = TRUE, control=list(maxit=2000))
MARSSparamCIs(kem, method = "parametric", nboot = 100)

dat <- unnamed[,2]
kem = MARSS(dat, model = list(U = matrix("U"), B = "unconstrained"), silent = TRUE, control=list(maxit=2000))
MARSSparamCIs(kem, method = "parametric", nboot = 100)


#figure 1 from 6-panel auto figure, do this to show pop estimates vs. observation counts
a <- crow
nyr = length(a[, 1])
plot(a[, 1], a[, 2], type = "p", bty = "L", xaxp = c(a[1, 
                                                       1], a[nyr, 1], nyr - 1), xlab = "", ylab = "Pop. Estimate", 
     ylim = c(0.9 * min(a[, 2], na.rm = TRUE), 1.1 * max(a[, 
                                                           2], na.rm = TRUE)), xlim = c(a[1, 1] - 1, a[nyr, 
                                                                                                       1] + 1))
lines(a[, 1], kem$states[1, ], col = "red")



#Try density-independent growth and PVA risk estiamtes for 14 stream segments

#doesn't really work at all, move on to next idea

#Model 3 subpopulations together, to find Q structure, B parameter


gaura <- read.csv("subcounts.csv", header = TRUE)
gaura.id <- as.factor(paste(gaura$Creek, gaura$Subunit, sep = "."))
pop.names <- levels(reorder(levels(gaura.id), c(8, 9, 10, 11, 12, 13, 14, 3, 4, 5, 6, 7, 1, 2)))

gaura.mat <- matrix(data = gaura$Count, nrow = 25)
colnames(gaura.mat) <- pop.names
rownames(gaura.mat) <- c(1989:2013)

gaura.dat <- log(t(gaura.mat + 1))

#First, test if structure between 3 subpopulations
unnamed <- log(rowSums(gaura.mat[, c(1:2)], na.rm = TRUE))
diamond <- log(rowSums(gaura.mat[, c(3:7)], na.rm = TRUE))
crow <- log(rowSums(gaura.mat[, c(7:14)], na.rm = TRUE))

pops <- cbind(unnamed, diamond, crow)

dat <- as.matrix(t(pops))


Z.model <- factor(c(1:3))
#three subpopulations makes more sense
Q.model <- c("unconstrained")
#unconstrained to see correlations of process error
U.model <- c("unequal")
#makes sense for trend to be unequal
R.model <- "diagonal and equal"
#makes sense for R (obs error) to be equal
A.model <- "scaling"
#having lower A for herbivory years better
B.models <- c("identity", "diagonal and equal", "diagonal and unequal")
#makes sense to be diagonal and unequal (density dependence)
x0.model <- "unequal"
V0.model <- "zero"
model.constant=list(U=U.model, Z=Z.model, R=R.model, A=A.model, Q=Q.model, x0=x0.model, V0=V0.model, tinitx=1)

#data mining to try it out
out.tab=NULL
fits=list()
for(B.model in B.models){
  fit.model = c(list(B=B.model), model.constant)
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


plot(1989:2013, fits[[2]]$states[1,], type = "l", col = "red", xlim = c(1989, 2013), ylim = c(3, 9), xlab = "Year", ylab = "Population")
lines(1989:2013, fits[[2]]$states[2,], col = "blue")
lines(1989:2013, fits[[2]]$states[3,], col = "green")
points(1989:2013, dat[1,], col = "red")
points(1989:2013, dat[2,], col = "blue")
points(1989:2013, dat[3,], col = "green")
lines(1989:2013, dat[1,], col = "red", lty = 2)
lines(1989:2013, dat[2,], col = "blue", lty = 2)
lines(1989:2013, dat[3,], col = "green", lty = 2)




######
#Now try different versions of trend, with 2000-2006 different for drought
#Also, try A with 2007-08 as greater observation error
A.time <- matrix(c("r","r","r"),3,1)
herb <- matrix(c("h", "h", "h"),3,1)
At <- array(A.time, dim = c(dim(A.time), dim(dat)[2]))
At[,,19:20] <- herb
A <- matrix(list(0, 0, 0),3,1)
A.models <- list(At, A)
names(A.models) <- c("time", "constant")


U.time <- matrix(c("u1", "u2", "u3"), 3, 1)
drought <- matrix(c("d1", "d2", "d3"), 3,1)
Ut <- array(U.time, dim = c(dim(U.time), dim(dat)[2]))
Ut[,,12:18] <- drought 
U <- matrix(list("u1","u2","u3"),3,1)
U.models <- list(Ut,U)
names(U.models) <- c("time", "constant")


Z.model <- factor(c(1:3))
#three subpopulations makes more sense
Q.model <- c("unconstrained")
#unconstrained to see correlations of process error
R.model <- "diagonal and equal"
#makes sense for R (obs error) to be equal
#having lower A for herbivory years better
B.model <- c("diagonal and equal")
#makes sense to be diagonal and unequal (density dependence)
x0.model <- "unequal"
V0.model <- "zero"
model.constant=list(Z=Z.model, R=R.model, B=B.model, Q=Q.model, x0=x0.model, V0=V0.model, tinitx=1)

#data mining to try it out
out.tab=NULL
fits=list()
for(i in 1:2){
  for(j in 1:2){
  fit.model = c(list(A=A.models[[i]], U = U.models[[j]]), model.constant)
  fit = MARSS(dat, model=fit.model,
              silent=TRUE, control=list(maxit=2000, safe = TRUE, trace = 1))
  out=data.frame(A=names(A.models)[i], U= names(U.models)[j],
                 logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                 num.iter=fit$numIter, converged=!fit$convergence)
  out.tab=rbind(out.tab,out)
  fits=c(fits,list(fit))
}
}

min.AICc=order(out.tab$AICc)
out.tab.1=out.tab[min.AICc,]
out.tab.1=cbind(out.tab.1, delta.AICc=out.tab.1$AICc-out.tab.1$AICc[1])
out.tab.1=cbind(out.tab.1, rel.like=exp(-1*out.tab.1$delta.AICc/2))
out.tab.1=cbind(out.tab.1, AIC.weight = out.tab.1$rel.like/sum(out.tab.1$rel.like))


plot(1989:2013, fits[[1]]$states[1,], type = "l", col = "red", xlim = c(1989, 2013), ylim = c(3, 9), xlab = "Year", ylab = "Population")
lines(1989:2013, fits[[1]]$states[2,], col = "blue")
lines(1989:2013, fits[[1]]$states[3,], col = "green")
points(1989:2013, dat[1,], col = "red")
points(1989:2013, dat[2,], col = "blue")
points(1989:2013, dat[3,], col = "green")
lines(1989:2013, dat[1,], col = "red", lty = 2)
lines(1989:2013, dat[2,], col = "blue", lty = 2)
lines(1989:2013, dat[3,], col = "green", lty = 2)


#########
#Should I try covariates?

load(file = "cov.pcpn.RData")
load(file = "cov.pdsi.RData")

covariates <- cov.pcpn
C.models <- c("unconstrained", "equal", "zero")
covs <- list(matrix(covariates[1,], nrow = 1), matrix(covariates[2,], nrow=1), matrix(covariates[3,], nrow=1), 
             matrix(covariates[4,], nrow=1))
######
#Also, try A with 2007-08 as greater observation error
A.time <- matrix(c("r","r","r"),3,1)
herb <- matrix(c("h", "h", "h"),3,1)
At <- array(A.time, dim = c(dim(A.time), dim(dat)[2]))
At[,,19:20] <- herb
A <- matrix(list(0, 0, 0),3,1)
A.models <- list(At, A)
names(A.models) <- c("time", "constant")

U.model <- c("unequal")
Z.model <- factor(c(1:3))
#three subpopulations makes more sense
Q.model <- c("unconstrained")
#unconstrained to see correlations of process error
R.model <- "diagonal and equal"
#makes sense for R (obs error) to be equal
#having lower A for herbivory years better
B.model <- c("diagonal and equal")
#makes sense to be diagonal and unequal (density dependence)
x0.model <- "unequal"
V0.model <- "zero"
model.constant=list(Z=Z.model, R=R.model, B=B.model, U=U.model, A = At, x0=x0.model, V0=V0.model, tinitx=1)

#data mining to try it out
out.tab=NULL
fits=list()

  for(C.model in C.models){
    for(j in 1:4){
    fit.model = c(list(C = C.model, c=covs[[j]]), model.constant)
    fit = MARSS(dat, model=fit.model,
                silent=TRUE, control=list(maxit=2000, safe = TRUE, trace = 1))
    out=data.frame(A=names(A.models)[i], C = C.model, c = j,
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   num.iter=fit$numIter, converged=!fit$convergence)
    out.tab=rbind(out.tab,out)
    fits=c(fits,list(fit))
  }
}


min.AICc=order(out.tab$AICc)
out.tab.1=out.tab[min.AICc,]
out.tab.1=cbind(out.tab.1, delta.AICc=out.tab.1$AICc-out.tab.1$AICc[1])
out.tab.1=cbind(out.tab.1, rel.like=exp(-1*out.tab.1$delta.AICc/2))
out.tab.1=cbind(out.tab.1, AIC.weight = out.tab.1$rel.like/sum(out.tab.1$rel.like))

#best precip model, but didnt' converge
fit.model = c(list(C="equal", Q = "diagonal and equal", c=covs[[1]]), model.constant)
fit <- MARSS(dat, model=fit.model,
             silent = TRUE, control=list(maxit = 5000, safe = TRUE, trace = 1))


plot(1989:2013, fit$states[1,], type = "l", col = "red", xlim = c(1989, 2013), ylim = c(3, 9), xlab = "Year", ylab = "Population")
lines(1989:2013, fit$states[2,], col = "blue")
lines(1989:2013, fit$states[3,], col = "green")
points(1989:2013, dat[1,], col = "red")
points(1989:2013, dat[2,], col = "blue")
points(1989:2013, dat[3,], col = "green")
lines(1989:2013, dat[1,], col = "red", lty = 2)
lines(1989:2013, dat[2,], col = "blue", lty = 2)
lines(1989:2013, dat[3,], col = "green", lty = 2)



#######################
#Time to try stream segments

gaura.mat <- matrix(data = gaura$Count, nrow = 25)
colnames(gaura.mat) <- pop.names
rownames(gaura.mat) <- c(1989:2013)

gaura.dat <- log(t(gaura.mat + 1))

######
#Also, try A with 2007-08 as greater observation error
A.time <- matrix(c("r","r","r"),14,1)
herb <- matrix(c("h", "h", "h"),14,1)
At <- array(A.time, dim = c(dim(A.time), dim(dat)[2]))
At[,,19:20] <- herb
A.model <- "zero"

U.model <-  matrix(list("uu","uu","ud","ud","ud","ud","ud","uc","uc","uc","uc","uc","uc","uc"),14,1)
Z.model <- factor(c(1:14))
#three subpopulations makes more sense
Q.model="unconstrained"
#unconstrained to see correlations of process error
R.model <- "diagonal and equal"
#makes sense for R (obs error) to be equal
#having lower A for herbivory years better
B1 <- c("identity")
B2=matrix(list(0),14,14)
diag(B2)=c("bu","bu","bd","bd","bd","bd","bd","bc","bc","bc","bc","bc","bc","bc")
#makes sense to be diagonal and unequal (density dependence)
x0.model <- "unequal"
V0.model <- "zero"
model.constant=list(Z=Z.model, Q=Q.model, R=R.model, B=B2, U=U.model, A = A.model, x0=x0.model, V0=V0.model, tinitx=1)

#data mining to try it out
out.tab=NULL
fits=list()

    fit.model = c(model.constant)
    fit = MARSS(gaura.dat, model=fit.model,
                silent=TRUE, control=list(maxit=2000, safe = TRUE, trace = 1))
    out=data.frame(
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   num.iter=fit$numIter, converged=!fit$convergence)
    out.tab=rbind(out.tab,out)
    fits=c(fits,list(fit))
  }
}


min.AICc=order(out.tab$AICc)
out.tab.1=out.tab[min.AICc,]
out.tab.1=cbind(out.tab.1, delta.AICc=out.tab.1$AICc-out.tab.1$AICc[1])
out.tab.1=cbind(out.tab.1, rel.like=exp(-1*out.tab.1$delta.AICc/2))
out.tab.1=cbind(out.tab.1, AIC.weight = out.tab.1$rel.like/sum(out.tab.1$rel.like))




#TRY SEGMENTS, but segment separately but with same covariate effects

covariates <- cov.pcpn
C.models <- c("equal", "zero")
covs <- list(matrix(covariates[1,], nrow = 1), matrix(covariates[2,], nrow=1), matrix(covariates[3,], nrow=1), 
             matrix(covariates[4,], nrow=1))
######
#Also, try A with 2007-08 as greater observation error
A.time <- matrix(rep("r", 14),14,1)
herb <- matrix(rep("h", 14),14,1)
At <- array(A.time, dim = c(dim(A.time), dim(dat)[2]))
At[,,19:20] <- herb
A.model <- "scaling"

U.model <- c("zero")
Z.model <- factor(c(1:14))
#three subpopulations makes more sense
Q.model <- c("diagonal and equal")
#unconstrained to see correlations of process error
R.model <- "diagonal and equal"
#makes sense for R (obs error) to be equal
#having lower A for herbivory years better
B.model <- c("identity")
#makes sense to be diagonal and unequal (density dependence)
x0.model <- "unequal"
V0.model <- "zero"
model.constant=list(Q=Q.model, Z=Z.model, R=R.model, B=B.model, U=U.model, A = A.model, x0=x0.model, V0=V0.model, tinitx=1)

#data mining to try it out
out.tab=NULL
fits=list()
for(C.model in C.models){
  for(j in 1:4){
    fit.model = c(list(C = C.model, c=covs[[j]]), model.constant)
    fit = MARSS(gaura.dat, model=fit.model,
                silent=TRUE, control=list(maxit=2000, safe = TRUE, trace = 1))
    out=data.frame(C = C.model, c = j,
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   num.iter=fit$numIter, converged=!fit$convergence)
    out.tab=rbind(out.tab,out)
    fits=c(fits,list(fit))
  }
}


min.AICc=order(out.tab$AICc)
out.tab.1=out.tab[min.AICc,]
out.tab.1=cbind(out.tab.1, delta.AICc=out.tab.1$AICc-out.tab.1$AICc[1])
out.tab.1=cbind(out.tab.1, rel.like=exp(-1*out.tab.1$delta.AICc/2))
out.tab.1=cbind(out.tab.1, AIC.weight = out.tab.1$rel.like/sum(out.tab.1$rel.like))



plot(1989:2013, fits[[1]]$states[1,], type = "l", col = "red", xlim = c(1989, 2013), ylim = c(0, 9), xlab = "Year", ylab = "Population")
lines(1989:2013, fits[[1]]$states[2,], col = "red")
points(1989:2013, gaura.dat[1,], col = "black")
points(1989:2013, gaura.dat[2,], col = "black")
lines(1989:2013, gaura.dat[1,], col = "black", lty = 2)
lines(1989:2013, gaura.dat[2,], col = "black", lty = 2)

plot(1989:2013, fits[[1]]$states[3,], type = "l", col = "red", xlim = c(1989, 2013), ylim = c(0, 9), xlab = "Year", ylab = "Population")
for(i in 4:7){
lines(1989:2013, fits[[1]]$states[i,], col = "red")
points(1989:2013, gaura.dat[i,], col = "black")
lines(1989:2013, gaura.dat[i,], col = "black", lty = 2)
}

plot(1989:2013, fits[[1]]$states[8,], type = "l", col = "red", xlim = c(1989, 2013), ylim = c(0, 7), xlab = "Year", ylab = "Population")
for(i in 9:14){
  lines(1989:2013, fits[[1]]$states[i,], col = "red")
  points(1989:2013, gaura.dat[i,], col = "black")
  lines(1989:2013, gaura.dat[i,], col = "black", lty = 2)
}


points(1989:2013, dat[1,], col = "red")
points(1989:2013, dat[2,], col = "blue")
points(1989:2013, dat[3,], col = "green")
lines(1989:2013, dat[2,], col = "blue", lty = 2)
lines(1989:2013, dat[3,], col = "green", lty = 2)

#########################
#One try at changing Z-matrix to make segments observations of 3 subpopulations


gaura.mat <- matrix(data = gaura$Count, nrow = 25)
colnames(gaura.mat) <- pop.names
rownames(gaura.mat) <- c(1989:2013)

gaura.dat <- log(t(gaura.mat + 1))

######
Z1 <- factor(c(rep("U",2), rep("D",5), rep("C",7)))
Z2 <- factor(c(1:14))

A.time <- matrix(rep("r", 14),14,1)
herb <- matrix(rep("h", 14),14,1)
At <- array(A.time, dim = c(dim(A.time), dim(dat)[2]))
At[,,19:20] <- herb
A.model <- "scaling"
U.model <-  "unequal"
Z.models <- list(Z1,Z2)
#three subpopulations makes more sense
Q.model=c("unconstrained")
#unconstrained to see correlations of process error
R.model <- "diagonal and equal"
#makes sense for R (obs error) to be equal
#having lower A for herbivory years better
B.model <- c("diagonal and unequal")
#makes sense to be diagonal and unequal (density dependence)
x0.model <- "unequal"
V0.model <- "zero"
model.constant=list(Q=Q.model, R=R.model, B=B.model, U=U.model, A = At, x0=x0.model, V0=V0.model, tinitx=1)

#data mining to try it out
out.tab=NULL
fits=list()
for(i in 1:2){
fit.model = c(list(Z=Z.models[[i]]), model.constant)
fit = MARSS(gaura.dat, model=fit.model,
            silent=TRUE, control=list(maxit=2000, safe = TRUE, trace = 1))
out=data.frame(Z=i,
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


unnamed <- log(rowSums(gaura.mat[, c(1:2)], na.rm = TRUE))
diamond <- log(rowSums(gaura.mat[, c(3:7)], na.rm = TRUE))
crow <- log(rowSums(gaura.mat[, c(7:14)], na.rm = TRUE))

pops <- cbind(unnamed, diamond, crow)


plot(1989:2013, fits[[1]]$states[1,], type = "l", col = "red", xlim = c(1989, 2013), ylim = c(2, 9), xlab = "Year", ylab = "Population")
lines(1989:2013, fits[[1]]$states[2,], col = "blue")
lines(1989:2013, fits[[1]]$states[3,], col = "green")
points(1989:2013, pops[,1], col = "red")
points(1989:2013, pops[,2], col = "blue")
points(1989:2013, pops[,3], col = "green")
lines(1989:2013, pops[,1], col = "red", lty = 2)
lines(1989:2013, pops[,2], col = "blue", lty = 2)
lines(1989:2013, pops[,3], col = "green", lty = 2)



######################
#What about extinction risk parameters now?
gaura.dat <- log(t(gaura.mat + 1))

#First, test if structure between 3 subpopulations
unnamed <- log(rowSums(gaura.mat[, c(1:2)], na.rm = TRUE))
diamond <- log(rowSums(gaura.mat[, c(3:7)], na.rm = TRUE))
crow <- log(rowSums(gaura.mat[, c(7:14)], na.rm = TRUE))

pops <- cbind(unnamed, diamond, crow)

dat <- as.matrix(t(pops))
######
#Now try different versions of trend, with 2000-2006 different for drought
#Also, try A with 2007-08 as greater observation error
A.time <- matrix(c("r","r","r"),3,1)
herb <- matrix(c("h", "h", "h"),3,1)
At <- array(A.time, dim = c(dim(A.time), dim(dat)[2]))
At[,,19:20] <- herb
A <- matrix(list(0, 0, 0),3,1)
A.models <- list(At, A)
names(A.models) <- c("time", "constant")

U.time <- matrix(c("u1", "u2", "u3"), 3, 1)
drought <- matrix(c("d1", "d2", "d3"), 3,1)
Ut <- array(U.time, dim = c(dim(U.time), dim(dat)[2]))
Ut[,,12:18] <- drought 
U <- matrix(list("u1","u2","u3"),3,1)
U.models <- list(Ut,U)
names(U.models) <- c("time", "constant")

Z.model <- factor(c(1:3))
#three subpopulations makes more sense
Q.model <- c("unconstrained")
#unconstrained to see correlations of process error
R.model <- "identity"
#makes sense for R (obs error) to be equal
#having lower A for herbivory years better
B.model <- c("identity")
x0.model <- "unequal"
V0.model <- "zero"
model.constant=list(Z=Z.model, R=R.model, B=B.model, Q=Q.model, x0=x0.model, V0=V0.model, tinitx=1)

#data mining to try it out
out.tab=NULL
fits=list()
# for(i in 1:2){
  for(j in 1:2){
    fit.model = c(list(A=A.models[[1]], U = U.models[[j]]), model.constant)
    fit = MARSS(dat, model=fit.model,
                silent=TRUE, control=list(maxit=2000, safe = TRUE, trace = 1))
    out=data.frame(A=names(A.models)[1], U= names(U.models)[j],
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   num.iter=fit$numIter, converged=!fit$convergence)
    out.tab=rbind(out.tab,out)
    fits=c(fits,list(fit))
  }
# }

min.AICc=order(out.tab$AICc)
out.tab.1=out.tab[min.AICc,]
out.tab.1=cbind(out.tab.1, delta.AICc=out.tab.1$AICc-out.tab.1$AICc[1])
out.tab.1=cbind(out.tab.1, rel.like=exp(-1*out.tab.1$delta.AICc/2))
out.tab.1=cbind(out.tab.1, AIC.weight = out.tab.1$rel.like/sum(out.tab.1$rel.like))


plot(1989:2013, fits[[1]]$states[1,], type = "l", col = "red", xlim = c(1989, 2013), ylim = c(3, 9), xlab = "Year", ylab = "Population")
lines(1989:2013, fits[[1]]$states[2,], col = "blue")
lines(1989:2013, fits[[1]]$states[3,], col = "green")
points(1989:2013, dat[1,], col = "red")
points(1989:2013, dat[2,], col = "blue")
points(1989:2013, dat[3,], col = "green")
lines(1989:2013, dat[1,], col = "red", lty = 2)
lines(1989:2013, dat[2,], col = "blue", lty = 2)
lines(1989:2013, dat[3,], col = "green", lty = 2)




#Herbivory accounted for or not, drought or not.
#extinction risk plots with basic calculations, just mu and sig2
#from powerpoint slide 8
U.1       0.0409
U.2       0.0111
U.3      -0.0733
Q.(1,1)   0.0978
Q.(2,1)   0.0973
Q.(3,1)   0.1181
Q.(2,2)   0.0971
Q.(3,2)   0.1177
Q.(3,3)   0.1430
x0.1      6.8153
x0.2      7.9021
x0.3      7.2654




CSEGtmufigure(N = 25, u = 0.0409, s2p = 0.0978)
CSEGtmufigure(N = 25, u = 0.0111, s2p = 0.0971)
CSEGtmufigure(N = 25, u = -0.0733, s2p = 0.1430)







