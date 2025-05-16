# clear environment
rm(list=ls())
# load libraries and utility codes
source("code/gen_data.R")
library(Rsurrogate)
library(TwoPhaseReg)

########################################
### Sim settings ######
########################################

### Sims settings
#m1 is MCAR
#m2 is MAR
#m3 is MNAR
setting.set = 1
# seed
sim_seed = 4
# set sample size
n1 = 1000 ### control group size
n0 = 1000 ### treatment group size



# set seed
set.seed(sim_seed)

# truth
truth = get.true(setting=1)


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
  Z = as.factor(rep(c(0,1), times = c(n0, n1))),
  Smiss = c(s0,s1)
)

# Simulate non-missingness indicators ###################
m2.1 = rep(1,n1)
m2.1[y1 < median(y1)] = rbinom(sum(y1 < median(y1)),1,0.70)
m2.0 = rep(1,n0)
m2.0[y0 > median(y0)] = rbinom(sum(y0 > median(y0)),1,0.60)
df$Miss = c(m2.0, m2.1)
df$Smiss[c(m2.0, m2.1) == 0] = NA

#known MAR weights
w.y1 = rep(1,n1)
w.y1[y1 < median(y1)] = 0.7
w.y0 = rep(1,n0)
w.y0[y0 > median(y0)] = 0.6
df$W = c(w.y0, w.y1)


y0.m = y0[m2.0 == 1]
s0.m = s0[m2.0 == 1]
w.y0.m = w.y0[m2.0 == 1]
y1.m = y1[m2.1 == 1]
s1.m = s1[m2.1 == 1]
w.y1.m = w.y1[m2.1 == 1]

##########################################################
#Estimates with complete data ############################
##########################################################

## Estimate R with parametric approach (Wang's approach) 
R_gold = Rsurrogate::R.s.estimate(sone=s1, 
                                  szero=s0, 
                                  yone = y1, 
                                  yzero = y0, 
                                  type = "model")

##########################################################
#Estimates with incomplete data ##########################
##########################################################

###########################################
## Layla's non-parametric
delta.ipw.nonparam = sum((1/w.y1[m2.1==1])*y1[m2.1==1])/length(y1)-sum((1/w.y0[m2.0==1])*y0[m2.0==1])/length(y0)
delta.s = delta.s.single(sone=s1[m2.1==1], szero=s0[m2.0==1], yone = y1[m2.1==1], yzero = y0[m2.0==1], weight.1 = w.y1[m2.1==1], weight.0 = w.y0[m2.0==1], n0.all = length(y0))

1 - delta.s/delta.ipw.nonparam

###########################################
## Estimate R with nonparametric approach (complete case)
# This will be biased, the estimator is overestimating the truth
R_miss_complete = Rsurrogate::R.s.estimate(sone= s1.m, 
                                           szero= s0.m, 
                                           yone = y1.m, 
                                           yzero =y0.m, 
                                           type = "model")

## from scratch
eq4 = lm(formula = y0.m ~ s0.m)
gamma0 = eq4$coefficients[1]
gamma1 = eq4$coefficients[2]
eq5 = lm(formula = y1.m ~ s1.m)
gamma2 = eq5$coefficients[1]
gamma3 = eq5$coefficients[2]
alpha0 = mean(s0.m) 
alpha1 = mean(s1.m)
delta = (gamma2 - gamma0)  + gamma3* alpha1 - + gamma1 *alpha0
delta
delta = mean(y1.m) - mean(y0.m)
delta
deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0
R = 1 - deltaS / delta
R

###########################################
## IPW
alpha0.w = sum((1/w.y0[m2.0==1])*s0[m2.0==1])/length(y0)
alpha1.w = sum((1/w.y1[m2.1==1])*s1[m2.1==1])/length(y1)
delta.ipw.nonparam = sum((1/w.y1[m2.1==1])*y1[m2.1==1])/length(y1)-sum((1/w.y0[m2.0==1])*y0[m2.0==1])/length(y0)
deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0.w 
R = 1 - deltaS / delta.ipw.nonparam
R

sim = data.frame(R = rep(NA,100),
                 delta = rep(NA,100),
                 delta.s = rep(NA,100) )
for(iter in 1:100){
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
    Z = as.factor(rep(c(0,1), times = c(n0, n1))),
    Smiss = c(s0,s1)
  )

  # Simulate non-missingness indicators ###################
  m2.1 = rep(1,n1)
  m2.1[y1 < median(y1)] = rbinom(sum(y1 < median(y1)),1,0.70)
  m2.0 = rep(1,n0)
  m2.0[y0 > median(y0)] = rbinom(sum(y0 > median(y0)),1,0.60)
  df$Miss = c(m2.0, m2.1)
  df$Smiss[c(m2.0, m2.1) == 0] = NA
  
  #known MAR weights
  w.y1 = rep(1,n1)
  w.y1[y1 < median(y1)] = 0.7
  w.y0 = rep(1,n0)
  w.y0[y0 > median(y0)] = 0.6
  df$W = c(w.y0, w.y1)
  
  
  y0.m = y0[m2.0 == 1]
  s0.m = s0[m2.0 == 1]
  w.y0.m = w.y0[m2.0 == 1]
  y1.m = y1[m2.1 == 1]
  s1.m = s1[m2.1 == 1]
  w.y1.m = w.y1[m2.1 == 1]
  eq4 = lm(formula = y0.m ~ s0.m, weights = 1/w.y0.m)
  gamma0 = eq4$coefficients[1]
  gamma1 = eq4$coefficients[2]
  eq5 = lm(formula = y1.m ~ s1.m,  weights = 1/w.y1.m)
  gamma2 = eq5$coefficients[1]
  gamma3 = eq5$coefficients[2]
  lm_alpha = lm(Smiss~Z, data = df, weights = 1/W)
  alpha0 = lm_alpha$coefficients[1]
  alpha1 = alpha0 + lm_alpha$coefficients[2]
  delta = (gamma2 - gamma0)  + gamma3* alpha1 - + gamma1 *alpha0
  deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0
  R = 1 - deltaS / delta
  sim$R[iter] = R
  sim$delta[iter] = delta
  sim$delta.s[iter] = deltaS
}
# work now!


###########################################
## SMLE
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


lm_alpha = lm(Smiss~Z, data = df, weights = 1/weights)
alpha0 = lm_alpha$coefficients[1]
alpha1 = alpha0 + lm_alpha$coefficients[2]
#### Calculate PTE 
delta = (gamma2 - gamma0)  + gamma3* alpha1 - + gamma1 *alpha0
deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0
1 - deltaS / delta
delta.ipw.nonparam = sum((1/w.y1[m2.1==1])*y1[m2.1==1])/length(y1)-sum((1/w.y0[m2.0==1])*y0[m2.0==1])/length(y0)
1 - deltaS / delta.ipw.nonparam


sim = data.frame(R = rep(NA,100),
                 delta = rep(NA,100),
                 delta.s = rep(NA,100),
                 delta.s.complete = rep(NA,100),
                 delta.y =  rep(NA, 100),
                 delta.complete =  rep(NA, 100),
                 R.y = rep(NA, 100), 
                 R.complete = rep(NA, 100))

for(iter in 1:100){
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
  
  # Simulate non-missingness indicators ###################
  m2.1 = rep(1,n1)
  m2.1[y1 < median(y1)] = rbinom(sum(y1 < median(y1)),1,0.70)
  m2.0 = rep(1,n0)
  m2.0[y0 > median(y0)] = rbinom(sum(y0 > median(y0)),1,0.60)
  df$Miss = c(m2.0, m2.1)
  df$Smiss[c(m2.0, m2.1) == 0] = NA
  
  #known MAR weights
  w.y1 = rep(1,n1)
  w.y1[y1 < median(y1)] = 0.7
  w.y0 = rep(1,n0)
  w.y0[y0 > median(y0)] = 0.6
  df$W = c(w.y0, w.y1)
  
  
  y0.m = y0[m2.0 == 1]
  s0.m = s0[m2.0 == 1]
  w.y0.m = w.y0[m2.0 == 1]
  y1.m = y1[m2.1 == 1]
  s1.m = s1[m2.1 == 1]
  w.y1.m = w.y1[m2.1 == 1]
  
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
  
  
  lm_alpha = lm(Smiss~Z, data = df, weights = 1/W)
  alpha0 = lm_alpha$coefficients[1]
  alpha1 = alpha0 + lm_alpha$coefficients[2]
  #### Calculate PTE 
  delta = (gamma2 - gamma0)  + gamma3* alpha1 - + gamma1 *alpha0
  deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0
  sim$delta[iter] = delta
  sim$delta.s[iter] = deltaS
  sim$R[iter] = 1 - deltaS / delta
  delta.ipw.nonparam = sum((1/w.y1[m2.1==1])*y1[m2.1==1])/length(y1)-sum((1/w.y0[m2.0==1])*y0[m2.0==1])/length(y0)
  sim$delta.y[iter] = delta.ipw.nonparam
  sim$R.y[iter] = 1 - deltaS / delta.ipw.nonparam
  alpha0 = mean(s0.m)
  alpha1 = mean(s1.m)
  delta = (gamma2 - gamma0)  + gamma3* alpha1 - + gamma1 *alpha0
  deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0
  sim$delta.complete[iter] = delta
  sim$delta.s.complete[iter] = deltaS
  sim$R.complete[iter] = 1 - deltaS / delta
}

# It works but not as good as IPW 







