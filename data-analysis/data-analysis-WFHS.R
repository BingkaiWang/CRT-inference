rm(list = ls())
set.seed(123)
library(tidyverse)
library(lme4)
library(geepack) #GEE
library(CRTgeeDR) #Aug-GEE
library(SuperLearner) # TMLE
library(tmle) # TMLE
SLmethods <- c("SL.glm", "SL.rpart", "SL.nnet")

# data preprocessing -----
d <- read.csv("data-analysis/36158-0001-Data.csv") %>%
  dplyr::select(ADMINLINK, WAVE, EMPLOYEE, STUDYGROUP, CONDITION, RMZFN, RMZEMP, SCWM_CWH) %>%
  mutate(SCWM_CWH = ifelse(SCWM_CWH >=0, SCWM_CWH, NA)) %>%
  filter(EMPLOYEE == 1 & WAVE %in% c(1,2)) %>%
  mutate(WAVE = ifelse(WAVE==1, "baseline", "Y")) %>%
  pivot_wider(names_from = WAVE, values_from = SCWM_CWH) %>%
  mutate(treatment = ifelse(CONDITION==1, 1, 0), cluster_id = as.factor(STUDYGROUP)) %>%
  dplyr::select(A=treatment, cluster_id, RMZFN, RMZEMP, baseline, Y) %>%
  filter(!is.na(Y))
d <- d %>% left_join(group_by(d, cluster_id) %>% summarise(N = n()), by = "cluster_id") 
d$A <- as.integer(d$A)
d$baseline[is.na(d$baseline)] <- mean(d$baseline, na.rm = T)
d$RMZFN[is.na(d$RMZFN)] <- median(d$RMZFN, na.rm = T)
d$RMZEMP[is.na(d$RMZEMP)] <- mean(d$RMZEMP, na.rm = T)
summary(d)
observed_data <- d
cl_data <- observed_data %>% group_by(cluster_id) %>% summarise_all(mean)
m <- length(unique(cl_data$cluster_id))
pi <- 0.5

# # cluster-level effects -------
# Unadjusted-Bugni
unadj_est <- mean(cl_data$Y[cl_data$A==1]) - mean(cl_data$Y[cl_data$A==0])
unadj_se <- sqrt(var(cl_data$Y[cl_data$A==1])/pi + var(cl_data$Y[cl_data$A==0])/(1-pi))/sqrt(m-1)

# gee
gee <- geeglm(Y~.-cluster_id, id = cluster_id, 
              data =observed_data, family = gaussian, corstr = "exchangeable")
gee_est <- summary(gee)$coefficients['A', 'Estimate']
gee_se <- summary(gee)$coefficients['A', 'Std.err'] * sqrt(m)/sqrt(m-4)

# glmm
glmm <- lmer(Y~(1|cluster_id) + .-cluster_id, data = observed_data)
glmm_est <- summary(glmm)$coefficients['A', 'Estimate']
glmm_se <- summary(glmm)$coefficients['A', 'Std. Error'] * sqrt(m)/sqrt(m-4)

# Aug-GEE
aug_gee <- geeDREstimation(Y~ A, id = "cluster_id", data = as.data.frame(observed_data), family = "gaussian",corstr = "exchangeable",
                           nameY = "Y", nameTRT = "A", model.augmentation.trt = Y ~ N+RMZFN+RMZEMP+baseline,
                           model.augmentation.ctrl = Y ~ N+RMZFN+RMZEMP+baseline)
aug_gee_est <- summary(aug_gee)$beta[2]
aug_gee_se <- summary(aug_gee)$se.robust[2] * sqrt(m)/sqrt(m-4)

# TMLE
set.seed(123)
cl_tmle <- tmle(Y = cl_data$Y, A = cl_data$A, W = cl_data[,-c(1,2,6)],
                Q.SL.library = SLmethods, g1W = rep(pi, m), family = "gaussian")
cl_est <- cl_tmle$estimates$ATE$psi
cl_se <- sqrt(cl_tmle$estimates$ATE$var.psi) * sqrt(m)/sqrt(m-4)

# SAAE-PM
kappa1 <- glm(A~N, family = gaussian, data = cl_data) %>% predict(type = "response")
kappa0 <- 1 - kappa1
zeta1 <- predict(lm(Y~N, data = cl_data[cl_data$A==1,]), newdata = cl_data, type = "response")
zeta0 <- predict(lm(Y~N, data = cl_data[cl_data$A==0,]), newdata = cl_data, type = "response")
eta <- observed_data %>%
  mutate(eta1= predict(glm(Y~., data = dplyr::select(observed_data[observed_data$A==1,],-c('A', 'cluster_id')), family = gaussian), newdata = observed_data, type = "response")) %>%
  mutate(eta0= predict(glm(Y~., data = dplyr::select(observed_data[observed_data$A==0,],-c('A', 'cluster_id')), family = gaussian), newdata = observed_data, type = "response")) %>%
  dplyr::select(cluster_id, eta1, eta0) %>%
  group_by(cluster_id) %>%
  summarise_all(mean) %>%
  as.data.frame()
mu_C1 <- cl_data$A / pi * (cl_data$Y - eta$eta1) + kappa1 / pi * (eta$eta1 - zeta1) + zeta1
mu_C0 <- (1-cl_data$A) / (1-pi) * (cl_data$Y - eta$eta0) + kappa0 / (1-pi) * (eta$eta0 - zeta0) + zeta0
saae_pwm_est_C <- mean(mu_C1) - mean(mu_C0)
saae_pwm_se_C <- sd(mu_C1 - mu_C0)/sqrt(m-4)

# # SAAE-ML
set.seed(123)
kappa1 <- SuperLearner(cl_data$A, X = cl_data[,c("N")], family = "gaussian",
                       SL.library = SLmethods, cvControl = list(V=5))$SL.predict
kappa0 <- 1 - kappa1
zeta_fit1 <- SuperLearner(cl_data$Y[cl_data$A==1], X = cl_data[cl_data$A==1,c("N")], family = "gaussian", 
                          SL.library = SLmethods, cvControl = list(V=2))
zeta_fit0 <- SuperLearner(cl_data$Y[cl_data$A==0], X = cl_data[cl_data$A==0,c("N")], family = "gaussian", 
                          SL.library = SLmethods, cvControl = list(V=2))
zeta1 <- predict(zeta_fit1, newdata = cl_data[,c("N")])$pred
zeta0 <- predict(zeta_fit0, newdata = cl_data[,c("N")])$pred
eta_fit1 <- SuperLearner(observed_data$Y[observed_data$A==1], X = dplyr::select(observed_data[observed_data$A==1,], -c("Y", "cluster_id","A")), family = "gaussian",
                         id = observed_data$cluster_id[observed_data$A==1], SL.library = SLmethods, cvControl = list(V=2))
eta_fit0 <- SuperLearner(observed_data$Y[observed_data$A==0], X = dplyr::select(observed_data[observed_data$A==0,], -c("Y", "cluster_id","A")), family = "gaussian",
                         id = observed_data$cluster_id[observed_data$A==0], SL.library = SLmethods, cvControl = list(V=2))
eta <- observed_data %>%
  mutate(eta1= predict(eta_fit1, newdata = dplyr::select(observed_data, -c("Y", "cluster_id","A")))$pred) %>%
  mutate(eta0 = predict(eta_fit0, newdata = dplyr::select(observed_data, -c("Y", "cluster_id",'A')))$pred) %>%
  dplyr::select(cluster_id, eta1, eta0) %>%
  group_by(cluster_id) %>%
  summarise_all(mean) %>%
  as.data.frame()
mu_C1 <- cl_data$A / pi * (cl_data$Y - eta$eta1) + kappa1 / pi * (eta$eta1 - zeta1) + zeta1
mu_C0 <- (1-cl_data$A) / (1-pi) * (cl_data$Y - eta$eta0) + kappa0 / (1-pi) * (eta$eta0 - zeta0) + zeta0
saae_ml_est_C <- mean(mu_C1) - mean(mu_C0)
saae_ml_se_C <- sd(mu_C1 - mu_C0)/sqrt(m-4)

summary_results <- matrix(data = c(unadj_est, gee_est, glmm_est, aug_gee_est, cl_est, saae_pwm_est_C, saae_ml_est_C,
                                   unadj_se, gee_se, glmm_se, aug_gee_se, cl_se, saae_pwm_se_C, saae_ml_se_C), nrow = 7, ncol = 2,
                          dimnames = list(c("unadj", "GEE", "LMM", "Aug-GEE", "TMLE", "PM","ML"),
                                          c("est", "sd"))) %>% as.data.frame %>%
  mutate(ci.lower = est + sd * qt(0.025, m), ci.upper = est + sd * qt(0.975, m), rov = sd^2/sd[1]^2)
round(summary_results,3)
xtable::xtable(summary_results[,-2])


# individual-level effect ----
# unadj
unadj_mu1 <- sum(cl_data$Y * cl_data$A * cl_data$N)/sum(cl_data$A * cl_data$N)
unadj_mu0 <- sum(cl_data$Y * (1-cl_data$A) * cl_data$N)/sum((1-cl_data$A) * cl_data$N)
unadj_est <- unadj_mu1 - unadj_mu0
unadj_se <- sd(cl_data$N/mean(cl_data$N) * (cl_data$A*(cl_data$Y-unadj_mu1)/pi - (1-cl_data$A)*(cl_data$Y - unadj_mu0)/(1-pi)))/sqrt(m-1)

# gee
gee <- geeglm(Y~.-cluster_id, id = cluster_id, weights = observed_data$N,
              data =observed_data, family = gaussian, corstr = "exchangeable")
gee_est <- summary(gee)$coefficients['A', 'Estimate']
gee_se <- summary(gee)$coefficients['A', 'Std.err'] * sqrt(m)/sqrt(m-4)

# glmm
scaled_data <- map_dfr(1:m, function(j){
  dd <- observed_data[observed_data$cluster_id == j,]
  N_j <- dd$N[1]
  map_dfr(1: N_j, function(k){
    mutate(dd, cluster_id = as.numeric(paste0(j, k)))
  })
})
glmm <- lmer(Y~(1|cluster_id) + . - cluster_id, REML=F, data = scaled_data)
glmm_est <- summary(glmm)$coefficients['A', 'Estimate']
glmm_se <- summary(glmm)$coefficients['A', 'Std. Error'] * sqrt(nrow(scaled_data)) / sqrt(nrow(observed_data)) * sqrt(m)/sqrt(m-4)

# aug-GEE
aug_gee <- geeDREstimation(Y~ A, id = "cluster_id", data = as.data.frame(scaled_data), family = "gaussian",corstr = "exchangeable",
                           nameY = "Y", nameTRT = "A", Y ~ N+RMZFN+RMZEMP+baseline,
                           model.augmentation.ctrl = Y ~ N+RMZFN+RMZEMP+baseline)
aug_gee_est <- summary(aug_gee)$beta[2]
aug_gee_se <- summary(aug_gee)$se.robust[2] * sqrt(nrow(scaled_data)) / sqrt(nrow(observed_data)) * sqrt(m)/sqrt(m-4)

#TMLE
set.seed(123)
I_tmle <- tmle(Y = observed_data$Y, A = observed_data$A, W = observed_data[,-c(1,2,6)], 
               Q.SL.library = SLmethods, g1W = rep(pi, nrow(observed_data)), family = "gaussian")
tmle_est <- I_tmle$estimates$ATE$psi
tmle_se <- cbind(IC = I_tmle$estimates$IC$IC.ATE, N = observed_data$N, id = observed_data$cluster_id) %>%
  as.data.frame %>% group_by(id) %>% summarise_all(sum) %>% dplyr::select(IC) %>% var %>% sqrt(.)/ nrow(observed_data) * sqrt(m-4)

# SAAE-PM
kappa1 <- glm(A~N, family = gaussian, data = cl_data) %>% predict(type = "response")
kappa0 <- 1 - kappa1
zeta1 <- predict(lm(Y~N, data = cl_data[cl_data$A==1,]), newdata = cl_data, type = "response")
zeta0 <- predict(lm(Y~N, data = cl_data[cl_data$A==0,]), newdata = cl_data, type = "response")
eta <- observed_data %>%
  mutate(eta1= predict(glm(Y~., data = dplyr::select(observed_data[observed_data$A==1,],-c('A', 'cluster_id')), family = gaussian), newdata = observed_data, type = "response")) %>%
  mutate(eta0= predict(glm(Y~., data = dplyr::select(observed_data[observed_data$A==0,],-c('A', 'cluster_id')), family = gaussian), newdata = observed_data, type = "response")) %>%
  dplyr::select(cluster_id, eta1, eta0) %>%
  group_by(cluster_id) %>%
  summarise_all(mean) %>%
  as.data.frame()
mu_C1 <- cl_data$A / pi * (cl_data$Y - eta$eta1) + kappa1 / pi * (eta$eta1 - zeta1) + zeta1
mu_C0 <- (1-cl_data$A) / (1-pi) * (cl_data$Y - eta$eta0) + kappa0 / (1-pi) * (eta$eta0 - zeta0) + zeta0
saae_pwm_est_I <- (mean(cl_data$N * mu_C1) - mean(cl_data$N * mu_C0))/mean(cl_data$N)
saae_pwm_se_I <- sqrt(var(cl_data$N *(mu_C1 - mu_C0 - saae_pwm_est_I)/mean(cl_data$N))/(m-4))

# # SAAE-ML
set.seed(123)
kappa1 <- SuperLearner(cl_data$A, X = cl_data[,c("N")], family = "gaussian",
                       SL.library = SLmethods, cvControl = list(V=5))$SL.predict
kappa0 <- 1 - kappa1
zeta_fit1 <- SuperLearner(cl_data$Y[cl_data$A==1], X = cl_data[cl_data$A==1,c("N")], family = "gaussian", 
                          SL.library = SLmethods, cvControl = list(V=2))
zeta_fit0 <- SuperLearner(cl_data$Y[cl_data$A==0], X = cl_data[cl_data$A==0,c("N")], family = "gaussian", 
                          SL.library = SLmethods, cvControl = list(V=2))
zeta1 <- predict(zeta_fit1, newdata = cl_data[,c("N")])$pred
zeta0 <- predict(zeta_fit0, newdata = cl_data[,c("N")])$pred
eta_fit1 <- SuperLearner(observed_data$Y[observed_data$A==1], X = dplyr::select(observed_data[observed_data$A==1,], -c("Y", "cluster_id","A")), family = "gaussian",
                         id = observed_data$cluster_id[observed_data$A==1], SL.library = SLmethods, cvControl = list(V=2))
eta_fit0 <- SuperLearner(observed_data$Y[observed_data$A==0], X = dplyr::select(observed_data[observed_data$A==0,], -c("Y", "cluster_id","A")), family = "gaussian",
                         id = observed_data$cluster_id[observed_data$A==0], SL.library = SLmethods, cvControl = list(V=2))
eta <- observed_data %>%
  mutate(eta1= predict(eta_fit1, newdata = dplyr::select(observed_data, -c("Y", "cluster_id","A")))$pred) %>%
  mutate(eta0 = predict(eta_fit0, newdata = dplyr::select(observed_data, -c("Y", "cluster_id",'A')))$pred) %>%
  dplyr::select(cluster_id, eta1, eta0) %>%
  group_by(cluster_id) %>%
  summarise_all(mean) %>%
  as.data.frame()
mu_C1 <- cl_data$A / pi * (cl_data$Y - eta$eta1) + kappa1 / pi * (eta$eta1 - zeta1) + zeta1
mu_C0 <- (1-cl_data$A) / (1-pi) * (cl_data$Y - eta$eta0) + kappa0 / (1-pi) * (eta$eta0 - zeta0) + zeta0
saae_ml_est_I <- (mean(cl_data$N * mu_C1) - mean(cl_data$N * mu_C0))/mean(cl_data$N)
saae_ml_se_I <- sqrt(var(cl_data$N *(mu_C1 - mu_C0 - saae_ml_est_I)/mean(cl_data$N))/(m-4))

summary_results <- matrix(data = c(unadj_est, gee_est, glmm_est, aug_gee_est, tmle_est, saae_pwm_est_I, saae_ml_est_I,
                                   unadj_se, gee_se, glmm_se, aug_gee_se, tmle_se, saae_pwm_se_I, saae_ml_se_I), nrow = 7, ncol = 2,
                          dimnames = list(c("unadj", "GEE", "LMM", "Aug-GEE", "TMLE", "PM","ML"),
                                          c("est", "sd"))) %>% as.data.frame %>%
  mutate(ci.lower = est + sd * qt(0.025, m), ci.upper = est + sd * qt(0.975, m), rov = sd^2/sd[1]^2)
round(summary_results,2)
xtable::xtable(summary_results[,-2])
