# clear environment
rm(list=ls())
# load libraries and utility codes
source("~/Research/missing_surrogate/code/4. sim code/gen_data.R")
#devtools::install_github("dragontaoran/TwoPhaseReg")
library(TwoPhaseReg)
library(tictoc)
library(Rsurrogate)

########################################
### Sim settings ######
########################################

### Sims settings
#m1 is MCAR
#m2 is MAR
#m3 is MNAR
setting.set = 1
# seed
sim_seed = 30
# set num of reps
num_reps = 1000
# set sample size
n1 = 1000 ### control group size
n0 = 1000 ### treatment group size
# prop of missingness for each treatment group, from 0 to 1
prop_obs_s0 = 0.6
prop_obs_s1 = 0.7

# set seed
set.seed(sim_seed)

truth = get.true(setting=1)

### Store sims data
sim_result = data.frame(
  sim = paste(sim_seed, 1:num_reps, sep = "-"),
  n1 = n1, n0 = n0, prop_obs_s0 = prop_obs_s0, prop_obs_s1 = prop_obs_s1,
  R.true = truth$R.true, delta.true = truth$delta.true, deltat.s.true = truth$delta.s.true,
  R.gold = NA, R.complete = NA, R.ipw.ipw = NA, R.ipw.complete = NA,
  R.smle.ipw = NA, R.smle.complete = NA, R.layla = NA,
  delta.gold = NA, delta.complete = NA, delta.ipw.ipw = NA, delta.ipw.complete = NA,
  delta.smle.ipw = NA, delta.smle.complete = NA, delta.layla = NA,
  delta.s.gold = NA, delta.s.complete = NA, delta.s.ipw.ipw = NA, delta.s.ipw.complete = NA,
  delta.s.smle.ipw = NA, delta.s.smle.complete = NA, delta.s.layla = NA
)


# Generate data 
data = gen.data(setting=1, n1=n1, n0=n0) 
# Define vectors for outcomes/surrogates in untreated/treated 
s1 = data$s1
y1 = data$y1
s0 = data$s0
y0 = data$y0

### Reformat data into dataframe

df = data.frame(
  Y = c(y0, y1),
  S = c(s0,s1),
  Z = rep(c(0,1), times = c(n0, n1)),
  Smiss = c(s0,s1)
)

##########################################################
#Estimates with complete data ############################
##########################################################

## Estimate R with parametric approach (Wang's approach) 
R_gold = Rsurrogate::R.s.estimate(sone=s1, 
                                  szero=s0, 
                                  yone = y1, 
                                  yzero = y0, 
                                  type = "model")
sim_result$R.gold[iter] = R_gold$R.s
sim_result$delta.gold[iter] = R_gold$delta
sim_result$delta.s.gold[iter] = R_gold$delta.s

##########################################################
#Estimates with incomplete data ##########################
##########################################################


# Simulate non-missingness indicators ###################
m2.1 = rep(1,n1)
m2.1[y1 < median(y1)] = rbinom(sum(y1 < median(y1)),1,prop_obs_s1)
m2.0 = rep(1,n0)
m2.0[y0 > median(y0)] = rbinom(sum(y0 > median(y0)),1,prop_obs_s0)
df$Miss = c(m2.0, m2.1)
df$Smiss[c(m2.0, m2.1) == 0] = NA

# known MAR weights
w.y1 = rep(1,n1)
w.y1[y1 < median(y1)] = prop_obs_s1
w.y0 = rep(1,n0)
w.y0[y0 > median(y0)] = prop_obs_s0
df$W = c(w.y0, w.y1)

# unknown MAR weights
## visualising missingness
library(visdat)
library(broom)
vis_miss(df)

## find probability of observed: response propensity stratification
prop_score_lm <- glm(
                Miss ~ Z + Y,
                data = df,
                family = binomial()
              )
df$W_est <- predict(prop_score_lm, type = "response", data = df)
mean((df$W-df$W_est)^2)

library(splines)
### calibration model
calib_model <- lm(
  Smiss ~
    ns(Y, df = 4) * Z,
  # use log of `wait_minutes_actual_avg`
  data = df)
df_calib <- calib_model |>
  augment(
    data = df |>
      drop_na()
  )


### Estimate R with nonparametric approach (complete case) ##########################
R_miss_complete = Rsurrogate::R.s.estimate(sone=s1[m2.1==1], 
                                           szero=s0[m2.0==1], 
                                           yone = y1[m2.1==1], 
                                           yzero = y0[m2.0==1], 
                                           type = "model")
sim_result$R.complete[iter] = R_miss_complete$R.s
sim_result$delta.complete[iter] = R_miss_complete$delta # this is the same as difference in mean outcomes
sim_result$delta.s.complete[iter] = R_miss_complete$delta.s

### IPW ##########################

### Estimate R with IPW only (semiparam)

#### Fit Y ~ Z + S
lm_ipw_semiparam = lm(formula = Y ~ Smiss*Z, 
                      data = df, 
                      weights = 1/W)
##### Get coefficients
gamma0 = as.numeric(lm_ipw_semiparam$coefficients[1])
gamma1 = as.numeric(lm_ipw_semiparam$coefficients[2])
gamma2 = gamma0 + as.numeric(lm_ipw_semiparam$coefficients[3])
gamma3 = gamma1 + as.numeric(lm_ipw_semiparam$coefficients[4])

#### Fit S ~ Z with IPW
lm_surr_weighted = lm(formula = Smiss ~ Z, data = df, 
                      weights = 1/W)
##### Get coefficients
alpha0_weighted = as.numeric(lm_surr_weighted$coefficients[1])
alpha1_weighted = alpha0_weighted + as.numeric(lm_surr_weighted$coefficients[2])

#### Calculate PTE 
delta = (gamma2 - gamma0) + gamma1 * (alpha1_weighted - alpha0_weighted) + (gamma3 - gamma1) * alpha1_weighted
deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0_weighted
sim_result$delta.ipw.ipw[iter] = delta
sim_result$delta.s.ipw.ipw[iter]  = deltaS
sim_result$R.ipw.ipw[iter]  = 1 - deltaS / delta 

#### Fit S ~ Z complete case
lm_surr_comp = lm(formula = Smiss ~ Z, data = df)
##### Get coefficients
alpha0_comp = as.numeric(coefficients(lm_surr_comp)[1])
alpha1_comp = alpha0_comp + as.numeric(coefficients(lm_surr_comp)[2])

#### Calculate PTE 
delta = (gamma2 - gamma0) + gamma1 * (alpha1_comp - alpha0_comp) + (gamma3 - gamma1) * alpha1_comp 
deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0_comp
sim_result$delta.ipw.complete[iter] = delta
sim_result$delta.s.ipw.complete[iter] = deltaS
sim_result$R.ipw.complete[iter]  = 1 - deltaS / delta 

## SMLE ################################################################################################

### Estimate R with smle & IPW (semiparam) ################################################

##### creating treatment group indicator
df$bs1 = as.numeric(df$Z == 0) ## Create bs1 = I(Z = 0), treated 
df$bs2 = as.numeric(df$Z == 1) ## Create bs2 = I(Z = 1), untreated

##### creating covariates for interaction model
df$Sinter = df$Smiss*df$Z

#### fit Y ~ Z*S with smle
fit = TwoPhaseReg::smle(Y = "Y", ## outcome  
                        X = c("Smiss", "Sinter"), ## "expensive" covariate (with missing data)
                        Z = "Z", ## "inexpensive" covariate (without missing data)
                        Bspline_Z = c("bs1", "bs2"), ## B-splines / indicators for Z 
                        data = df, ## name of the dataset 
                        noSE = FALSE, ## do you want standard error estimates? 
                        TOL = 1E-4, ## used to define convergence
                        MAX_ITER = 500, ## used to keep it from blowing up 
                        model = "linear" ## type of outcome model - we want OLS
)
##### Get coefficients
gamma0 = as.numeric(fit$coefficients[1])
gamma1 = as.numeric(fit$coefficients[2])
gamma2 = gamma0 + as.numeric(fit$coefficients[4])
gamma3 = gamma1 + as.numeric(fit$coefficients[3])

#### Calculate PTE 
delta = (gamma2 - gamma0) + gamma1 * (alpha1_weighted - alpha0_weighted) + (gamma3 - gamma1) * alpha1_weighted
deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0_weighted
sim_result$delta.smle.ipw[iter] = delta
sim_result$delta.s.smle.ipw[iter] = deltaS
sim_result$R.smle.ipw[iter]  = 1 - deltaS / delta  

#### Calculate PTE 
delta = (gamma2 - gamma0) + gamma1 * (alpha1_comp - alpha0_comp) + (gamma3 - gamma1) * alpha1_comp 
deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0_comp
sim_result$delta.smle.complete[iter] = delta
sim_result$delta.s.smle.complete[iter] = deltaS
sim_result$R.smle.complete[iter]  = 1 - deltaS / delta  
  










