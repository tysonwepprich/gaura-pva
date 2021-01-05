# COBP GLMM

# Setup----
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(lme4)
library(mgcv)

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
  filter(segment %!in% c("C-VII", "C-VIII")) %>%
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
  filter(segment %!in% c("C-VII", "C-VIII")) %>%
  spread(year, flower.count) %>% 
  dplyr::select(-segment) %>%
  data.matrix()

ydata[,14:30] <- polydat # some excel errors noted in segment sums compared to polygon totals


# covariates
temp <- read.csv("temp.cheyenne.csv", header = TRUE, skip = 1) %>% 
  dplyr::select(year, temp.summer.apr.aug) %>% 
  filter(year > 1986) %>% 
  mutate(ztemp = scale(temp.summer.apr.aug)) %>% 
  mutate(ztemp_lag1 = lag(ztemp, 1),
         ztemp_lag2 = lag(ztemp, 2))
prec <- read.csv("precip.cheyenne.csv", header = TRUE, skip = 1) %>% 
  dplyr::select(year, precip.summer.apr.aug) %>% 
  filter(year > 1986) %>% 
  mutate(zprec = scale(precip.summer.apr.aug)) %>% 
  mutate(zprec_lag1 = lag(zprec, 1),
         zprec_lag2 = lag(zprec, 2))
flow <- read.csv("flow_estimates.csv", header = TRUE) %>% 
  dplyr::select(year, flow) %>% 
  filter(year > 1986) %>% 
  mutate(zflow = log(flow)) %>% 
  mutate(zflow_lag1 = lag(zflow, 1),
         zflow_lag2 = lag(zflow, 2))


counts <- data.frame(ydata) %>% 
  mutate(segment = unique(segdat$segment)[-c(7:8)]) %>% 
  tidyr::gather(key = year, value = count, X1989:X2018) %>% 
  mutate(year = as.numeric(gsub(pattern = "X", replacement = "", x = year, fixed = TRUE))) %>% 
  left_join(temp) %>% 
  left_join(prec) %>% 
  left_join(flow)

# add log population size and log lambda
counts2 <- counts %>% 
  group_by(segment) %>% 
  arrange(year) %>% 
  mutate(count_lag1 = lag(count, 1),
         count_lag2 = lag(count, 2),
         logN = log(count + 1),
         logN_lag1 = log(count_lag1 + 1),
         logN_lag2 = log(count_lag2 + 1),
         # logN_lag1 = scale(log(count_lag1 + 1)),
         # logN_lag2 = scale(log(count_lag2 + 1)),
         loglam = logN - log(count_lag1 + 1)) %>% 
  ungroup() %>% 
  filter(complete.cases(.))


mod2 <- lmer(loglam ~
              # logN_lag1 +
              # logN_lag2 +
               count_lag1 +
               count_lag2 +
              ztemp +
              ztemp_lag1 +
              ztemp_lag2 +
              # zprec +
              # zprec_lag1 +
              # zprec_lag2 +
              zflow +
              zflow_lag1 +
              zflow_lag2 +
              # ztemp : zflow +
              # ztemp_lag1 : zflow_lag1 +
              # ztemp_lag2 : zflow_lag2 +
               (1|year),
            data = counts2)
summary(mod2)
MuMIn::r.squaredGLMM(mod2)

write.csv(broom::tidy(mod2), file = "glmm_results.csv")

# summer temp and precip strongly negatively correlated (-.7) so could just include one
library(corrplot)
M <- counts %>% 
  dplyr::select(ztemp:ztemp_lag2, zprec:zprec_lag2, zflow:zflow_lag2) %>% 
  as.matrix()
M<-cor(M)

corrplot(M, method="circle")
