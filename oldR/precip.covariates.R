#process precipitation data for MARSS covariates
library(plyr)
library(lubridate)
library(reshape)

setwd("C:/Users/Tyson/Dropbox/Gaura")
zndx <- read.csv("wy08zndx.csv", header = FALSE)
pdsi <- read.csv("wy08pdsi.csv", header = FALSE)
phdi <- read.csv("wy08phdi.csv", header = FALSE)
pmdi <- read.csv("wy08pmdi.csv", header = FALSE)
pcpn <- read.csv("wy08pcpn.csv", header = FALSE)

grow.month <- c("apr", "may", "jun", "jul", "aug", "sep")
flow.month <- c("jun", "jul", "aug")

pcpn$year <- colsplit(as.character(pcpn$V1), split="0801", c("test1", "test2"))[,2]
pcpn <- pcpn[,-1]
colnames(pcpn) <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec", "year")
pcpn.melt <- melt(pcpn, id = "year")
pcpn.season <- ddply(pcpn.melt, .(year), summarize,
                     pcpn.grow = sum(value[variable %in% grow.month]),
                     pcpn.flow = sum(value[variable %in% flow.month]))
pcpn.season$z.pcpn.flow = scale(pcpn.season$pcpn.flow, center = TRUE, scale = TRUE)
pcpn.season$z.pcpn.grow = scale(pcpn.season$pcpn.grow, center = TRUE, scale = TRUE)

pmdi$year <- colsplit(as.character(pmdi$V1), split="0808", c("test1", "test2"))[,2]
pmdi <- pmdi[,-1]
colnames(pmdi) <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec", "year")
pmdi.melt <- melt(pmdi, id = "year")
pmdi.season <- ddply(pmdi.melt, .(year), summarize,
                     pmdi.grow = sum(value[variable %in% grow.month]),
                     pmdi.flow = sum(value[variable %in% flow.month]))
pmdi.season$z.pmdi.flow = scale(pmdi.season$pmdi.flow, center = TRUE, scale = TRUE)
pmdi.season$z.pmdi.grow = scale(pmdi.season$pmdi.grow, center = TRUE, scale = TRUE)


phdi$year <- colsplit(as.character(phdi$V1), split="0806", c("test1", "test2"))[,2]
phdi <- phdi[,-1]
colnames(phdi) <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec", "year")
phdi.melt <- melt(phdi, id = "year")
phdi.season <- ddply(phdi.melt, .(year), summarize,
                     phdi.grow = sum(value[variable %in% grow.month]),
                     phdi.flow = sum(value[variable %in% flow.month]))
phdi.season$z.phdi.flow = scale(phdi.season$phdi.flow, center = TRUE, scale = TRUE)
phdi.season$z.phdi.grow = scale(phdi.season$phdi.grow, center = TRUE, scale = TRUE)


pdsi$year <- colsplit(as.character(pdsi$V1), split="0805", c("test1", "test2"))[,2]
pdsi <- pdsi[,-1]
colnames(pdsi) <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec", "year")
pdsi.melt <- melt(pdsi, id = "year")
pdsi.season <- ddply(pdsi.melt, .(year), summarize,
                     pdsi.grow = sum(value[variable %in% grow.month]),
                     pdsi.flow = sum(value[variable %in% flow.month]))
pdsi.season$z.pdsi.flow = scale(pdsi.season$pdsi.flow, center = TRUE, scale = TRUE)
pdsi.season$z.pdsi.grow = scale(pdsi.season$pdsi.grow, center = TRUE, scale = TRUE)


zndx$year <- colsplit(as.character(zndx$V1), split="0807", c("test1", "test2"))[,2]
zndx <- zndx[,-1]
colnames(zndx) <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec", "year")
zndx.melt <- melt(zndx, id = "year")
zndx.season <- ddply(zndx.melt, .(year), summarize,
                     zndx.grow = sum(value[variable %in% grow.month]),
                     zndx.flow = sum(value[variable %in% flow.month]))
zndx.season$z.zndx.flow = scale(zndx.season$zndx.flow, center = TRUE, scale = TRUE)
zndx.season$z.zndx.grow = scale(zndx.season$zndx.grow, center = TRUE, scale = TRUE)

pairs(cbind(pcpn.season[,4], pmdi.season[,4], phdi.season[,4], pdsi.season[,4], zndx.season[,4]))

#takeaway: precipitation most variable, weaker (.7) correlation with drought indices
#drought indices highly correlated (except for zndx)
#use either pdsi or pcpn as covariates

#precipitation with time lags
pcpn.lag0 <- pcpn.season$z.pcpn.grow[10:34]
pcpn.lag1 <- pcpn.season$z.pcpn.grow[9:33]
pcpn.lag2 <- pcpn.season$z.pcpn.grow[8:32]
pcpn.lag3 <- pcpn.season$z.pcpn.grow[7:31]
#pdsi
pdsi.lag0 <- pdsi.season$z.pdsi.grow[10:34]
pdsi.lag1 <- pdsi.season$z.pdsi.grow[9:33]
pdsi.lag2 <- pdsi.season$z.pdsi.grow[8:32]
pdsi.lag3 <- pdsi.season$z.pdsi.grow[7:31]

cov.pcpn  <- as.matrix(rbind(pcpn.lag0,pcpn.lag1, pcpn.lag2, pcpn.lag3))
cov.pdsi <- as.matrix(rbind(pdsi.lag0, pdsi.lag1, pdsi.lag2, pdsi.lag3))

save(cov.pcpn, file = "cov.pcpn.RData")
save(cov.pdsi, file = "cov.pdsi.RData")
