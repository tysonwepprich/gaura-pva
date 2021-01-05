# Environmental covariates

library(dplyr)
library(ggplot2)
library(mgcv)
library(lubridate)
theme_set(theme_bw(base_size = 14))
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

temp <- read.csv("temp.cheyenne.csv", skip = 1, header = TRUE)
prec <- read.csv("precip.cheyenne.csv", skip = 1, header = TRUE)
water <- read.csv("crow.flow.csv", skip = 1, header = TRUE)

waterdf <- data.frame(year = 1994:2018, temp = temp$temp.wateryear[60:84], 
                      prec = prec$precip.wateryear[60:84], flow = water$flow.wateryear[2:26])

mod <- gam(flow ~ ti(temp, k = 5) + ti(prec, k = 5) + ti(temp, prec), data = waterdf, family = nb(theta = NULL, link = "log"))
waterdf$pred <- predict(mod, type = "response")

ggplot(temp, aes(x = year, y = temp.wateryear)) + 
  geom_point() +
  geom_smooth(method = "gam")
ggplot(prec, aes(x = year, y = precip.wateryear)) + 
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))
ggplot(water, aes(x = year, y = flow.wateryear)) + 
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), method.args = list(family = "poisson"))

# pca?

pcdat <- log(water[2:25, 2:13])
pcdat <- temp[54:84, 2:13]
pcdat <- cbind(temp[54:84, 2:13], log(prec[54:84, 2:13] + 0.01))

pc <- prcomp(pcdat, center = TRUE, scale. = TRUE)
pc
summary(pc)


library(randomForest)

# predict monthly water, to then estimate annual wateryear flow

lags <- seq(18)
# lag_names <- paste("lag", formatC(lags, width = nchar(max(lags)), flag = "0"), 
                   # sep = "_")
# lag_functions <- setNames(paste("dplyr::lag(., ", lags, ")"), lag_names)

# d %>% mutate_at(vars(x), funs_(lag_functions))

lag_names <- paste("temp_lag", formatC(lags, width = nchar(max(lags)), flag = "0"), 
                   sep = "_")
lag_functions <- setNames(paste("dplyr::lag(., ", lags, ")"), lag_names)

temp_mon <- temp %>% 
  tidyr::gather(key = month, value = temp, jan:dec) %>% 
  mutate(date = lubridate::ymd(paste(year, month, 1))) %>% 
  dplyr::select(date, temp) %>% 
  arrange(date) %>% 
  mutate_at(vars(-date), funs_(lag_functions)) %>% 
  filter(complete.cases(.))

lag_names <- paste("prec_lag", formatC(lags, width = nchar(max(lags)), flag = "0"), 
                   sep = "_")
lag_functions <- setNames(paste("dplyr::lag(., ", lags, ")"), lag_names)

prec_mon <- prec %>% 
  tidyr::gather(key = month, value = prec, jan:dec) %>% 
  mutate(date = lubridate::ymd(paste(year, month, 1))) %>% 
  dplyr::select(date, prec) %>% 
  arrange(date) %>% 
  mutate(prec = log(prec + 0.1)) %>% 
  mutate_at(vars(-date), funs_(lag_functions)) %>% 
  filter(complete.cases(.))

water_mon <- water %>% 
  tidyr::gather(key = month, value = flow, jan:dec) %>% 
  mutate(date = lubridate::ymd(paste(year, month, 1))) %>% 
  dplyr::select(date, flow) %>% 
  arrange(date) %>% 
  filter(complete.cases(.)) %>% 
  left_join(temp_mon) %>% 
  left_join(prec_mon)

water_mon$month <- as.factor(month(water_mon$date))
water_mon$year <- year(water_mon$date)

mod <- randomForest(x = water_mon[, 3:42], y = log(water_mon[, 2]), ntree = 1000)
water_mon$pred <- exp(predict(mod, newdata = water_mon))


ggplot(water_mon, aes(x = pred, y = flow, color = month)) +
  geom_point() +
  facet_wrap(~year, scales = "free") +
  geom_abline(slope = 1, intercept = 0) 
# 
# # what about PCA first on vars?
# water_mon_pca <- prcomp(as.matrix(water_mon[3:40]), center = TRUE, scale. = TRUE)
# wpc <- as.data.frame(water_mon_pca$x[,1:20])
# wpc$month <- as.factor(month(water_mon$date))
# wpc$year <- year(water_mon$date)
# 
# mod2 <- randomForest(x = wpc, y = water_mon[, 2])
# water_mon$pred <- predict(mod2, newdata = water_mon)





newdata <- temp_mon %>% 
  left_join(prec_mon) %>% 
  mutate(month = as.factor(month(date)),
         year = year(date)) %>% 
  filter(year >= 1980) %>% 
  mutate(predflow = exp(predict(mod, newdata = .)))

pred_wateryr <- newdata %>% 
  mutate(pred_flow = c(rep(NA, 11), zoo::rollapply(predflow, 12, mean)),
         obs_prec = c(rep(NA, 11), zoo::rollapply(exp(prec) - .1, 12, sum)),
         obs_temp = c(rep(NA, 11), zoo::rollapply(temp, 12, mean))) %>% 
  filter(month == 9) %>% 
  dplyr::select(year:obs_temp) %>% 
  left_join(waterdf)

ggplot(data  = pred_wateryr, aes(x = year)) + 
  geom_point(aes(y = pred_flow), shape = 2) +
  # geom_point(aes(y = obs_prec), color = "red") +
  geom_point(aes(y = flow), color = "black") +
  xlab("Year") +
  ylab("Observed and predicted water-year mean flow")


# env cov with flow prediction only in years without observation
envcov <- pred_wateryr %>% 
  mutate(flow = ifelse(is.na(flow), pred_flow, flow)) %>% 
  dplyr::select(year, obs_prec, obs_temp, flow)

write.csv(envcov, "flow_estimates.csv")




