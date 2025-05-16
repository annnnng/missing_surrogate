# clear environment
rm(list=ls())
# load libraries and utility codes
source("code/gen_data.R")
#devtools::install_github("dragontaoran/TwoPhaseReg")
library(TwoPhaseReg)
library(tictoc)
#library(Rsurrogate)

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
num_reps=10
# set sample size
n1 = 1000 ### control group size
n0 = 1000 ### treatment group size
# prop of missingness for each treatment group, from 0 to 1
#prop_obs_y0 = 0.7
#prop_obs_y1 = 0.6
prop_obs_s0 = 0.7
prop_obs_s1 = 0.8

for (n1 in c(500, 1000, 1500)) {
  print(paste("System time =", Sys.time()))
  tic(paste("Sims with N1 =", n1)) 
  
  for (n0 in c(500, 1000, 1500)){
    
    # set seed
    set.seed(sim_seed)
    
    ### Store sims data
    sim_result = data.frame(
      sim = paste(sim_seed, 1:num_reps, sep = "-"),
      n1 = NA, n0 = NA, prop_obs_s0 = NA, prop_obs_s1 = NA,
      delta.s.true = NA, delta.true = NA, R.true = NA,
      betahat1 = NA, betahat1_star = NA, R.model = NA,
      betahat1_star.complete = NA, R.model.complete = NA,
      betahat1_star.ipw = NA, R.model.ipw = NA,
      betahat1_star.smle = NA, R.model.smle = NA
    )
    
    for (iter in 1:num_reps){
      ########################################
      ### Generate data: 
      # y - outcome, 
      # s - surrogate, 
      # z - treatment group (0 = treatment, 1 = control)
      data = gen.data(setting=setting.set, n1=n1, n0=n0)
      s1 = data$s1
      y1 = data$y1
      s0 = data$s0
      y0 = data$y0
      dat = data.frame(Y = c(y0, y1), 
                       Z = c(rep(0, n0), rep(1, n1)), 
                       S = c(s0, s1))
      
      sim_result[iter,c("n1","n0","prop_obs_s1","prop_obs_s0")] = c(n1,n0,prop_obs_s1,prop_obs_s0)
      
      ########################################
      ### Truth ######
      ########################################
      # return data gen settings, other setting will return empirical values
      sim_result[iter,c("delta.s.true","delta.true","R.true")] = get.true(setting=1) 
      
      ########################################
      #no missingness: robust estimate#
      ########################################
      # observed total treatment effect
      #delta.obs.robust = mean(y1)-mean(y0)
      #delta.s.robust = delta.s.single(sone=s1, szero=s0, yone = y1, yzero = y0)
      #R.robust = 1-delta.s.robust/delta.obs.robust
      
      # save to dataframe
      #sim_result$delta.s.robust[iter] = delta.s.robust
      #sim_result$delta.robust[iter] = delta.robust
      #sim_result$R.robust[iter] = R.robust
      
      ########################################
      #no missingness: model estimate#
      ########################################
      ## Fit Y ~ Z 
      modYgivZ = lm(formula = Y ~ Z, 
                    data = dat)
      betahat1 = coefficients(modYgivZ)[2]
      
      ## Fit Y ~ Z + S
      modYgivZS = lm(formula = Y ~ Z + S, 
                     data = dat)
      betahat1_star = coefficients(modYgivZS)[2]
      
      ## Calculate PTE = R_F = 1 - \beta_1^* / \beta_1
      R.model = 1 - (betahat1_star / betahat1)
      
      # save to dataframe
      sim_result$betahat1_star[iter] = betahat1_star
      sim_result$betahat1[iter] = betahat1
      sim_result$R.model[iter] = R.model
      
      ########################################
      ### Add missingness MAR ######
      ########################################
      # observed indicator for outcome: 1 = observed, 0 = missing
      #y1_obs_ind = rep(1,n1)
      #y1_obs_ind[y1 < median(y1)] = rbinom(sum(y1 < median(y1)),1,prop_obs_y1)
      #y0_obs_ind = rep(1,n0)
      #y0_obs_ind[y0 > median(y0)] = rbinom(sum(y0 > median(y0)),1,prop_obs_y0)
      # observed indicator for surrogate: 1 = observed, 0 = missing
      #s1_obs_ind = rep(1,n1)
      #s1_obs_ind[s1 < median(s1)] = rbinom(sum(s1 < median(s1)),1,prop_obs_s1)
      #s0_obs_ind = rep(1,n0)
      #s0_obs_ind[s0 > median(s0)] = rbinom(sum(s0 > median(s0)),1,prop_obs_s0)
      
      s1_obs_ind = rbinom(n1,1,0.70)
      s0_obs_ind = rbinom(n0,1,0.70)
      
      #known MAR weights
      #w.y1 = rep(1,n1)
      #w.y1[y1 < median(y1)] = prop_obs_y1
      #w.y0 = rep(1,n0)
      #w.y0[y0 > median(y0)] = prop_obs_y0
      #w.s1 = rep(1,n1)
      #w.s1[s1 < median(s1)] = prop_obs_s1
      #w.s0 = rep(1,n0)
      #w.s0[s0 > median(s0)] =i prop_obs_s0
      
      w.s1 = rep(1/0.7, n1)
      w.s0 = rep(1/0.7, n0)
      
      # create missing data
      obs_ind = c(s1_obs_ind, s0_obs_ind)
      dat$Smiss = dat$S
      dat$Smiss[obs_ind == 0] = NA
      
      
      ##########################################################
      # Robust - Complete case #
      ##########################################################
      #delta.obs.robust.complete = mean(y1[s1_obs_ind==1])-mean(y0[s0_obs_ind==1])
      #delta.s.robust.complete = delta.s.single(sone=s1[s1_obs_ind==1], szero=s0[s0_obs_ind==1], yone = y1[s1_obs_ind==1], yzero = y0[s0_obs_ind==1])
      #R.robust.complete = 1-delta.s.robust.complete/delta.obs.robust.complete
      
      # save to dataframe
      #sim_result$delta.obs.robust.complete[iter] = delta.obs.robust.complete
      #sim_result$delta.s.robust.complete[iter] = delta.s.robust.complete
      #sim_result$R.robust.complete[iter] = R.robust.complete
      
      ##########################################################
      # Robust - Inverse Proportional Weighted #
      ##########################################################
      #delta.obs.robust.weighted = sum((1/w.y1[y1_obs_ind==1])*y1[y1_obs_ind==1])/length(y1)-sum((1/w.y0[y0_obs_ind==1])*y0[y0_obs_ind==1])/length(y0)
      #delta.s.robust.weighted = delta.s.single(sone=s1[y1_obs_ind==1], szero=s0[y0_obs_ind==1], yone = y1[y1_obs_ind==1], yzero = y0[y0_obs_ind==1], weight.1 = w.y1[y1_obs_ind==1], weight.0 = w.y0[y0_obs_ind==1], n0.all = length(y0))
      #R.robust.weighted = 1-delta.s.robust.weighted/delta.obs.robust.weighted 
      
      # save to dataframe
      #sim_result$delta.obs.robust.weighted[iter] = delta.obs.robust.weighted
      #sim_result$delta.s.robust.weighted[iter] = delta.s.robust.weighted
      #sim_result$R.robust.weighted[iter] = R.robust.weighted
      
      ##########################################################
      # Parametric - Complete case #
      ##########################################################
      ## Fit Y ~ Z + S
      modYgivZS = lm(formula = Y ~ Z + Smiss, 
                     data = dat)
      betahat1_star.complete = coefficients(modYgivZS)[2]
      
      ## Calculate PTE = R_F = 1 - \beta_1^* / \beta_1
      R.model.complete = 1 - (betahat1_star.complete / betahat1)
      
      # save to dataframe
      sim_result$betahat1_star.complete[iter] = betahat1_star.complete
      sim_result$R.model.complete[iter] = R.model.complete
      
      
      ##########################################################
      # Parametric - Inverse Proportional Weighted  #
      ##########################################################
      #y_weights = c(1/w.y1[y1_obs_ind==1], 1/w.y0[y0_obs_ind==1])
      dat$s.W = NA
      dat$s.W[obs_ind == 1] = c(1/w.s1[s1_obs_ind==1], 1/w.s0[s0_obs_ind==1])
      
      ## Fit Y ~ Z + S
      modYgivZS = lm(formula = Y ~ Z + Smiss, 
                     data = dat,
                     weights = s.W)
      betahat1_star.ipw = coefficients(modYgivZS)[2]
      
      ## Calculate PTE = R_F = 1 - \beta_1^* / \beta_1
      R.model.ipw = 1 - (betahat1_star.ipw / betahat1)
      
      # save to dataframe
      sim_result$betahat1_star.ipw[iter] = betahat1_star.ipw
      sim_result$R.model.ipw[iter] = R.model.complete
      
      ##########################################################
      # Parametric - Ran's two-phase SMLE #
      ##########################################################
      # creating treatment group indicator
      dat$bs1 = as.numeric(dat$Z == 0) ## Create bs1 = I(Z = 0), treated 
      dat$bs2 = as.numeric(dat$Z == 1) ## Create bs2 = I(Z = 1), untreated
      
      fit = TwoPhaseReg::smle(Y = "Y", ## outcome  
                              X = "Smiss", ## "expensive" covariate (with missing data)
                              Z = "Z", ## "inexpensive" covariate (without missing data)
                              Bspline_Z = c("bs1", "bs2"), ## B-splines / indicators for Z 
                              data = dat, ## name of the dataset 
                              noSE = FALSE, ## do you want standard error estimates? 
                              TOL = 1E-4, ## used to define convergence
                              MAX_ITER = 500, ## used to keep it from blowing up 
                              model = "linear" ## type of outcome model - we want OLS
      )
      betahat1_star.smle = fit$coefficients[3]
      
      ## Calculate PTE = R_F = 1 - \beta_1^*(miss) / \beta_1
      R.model.smle = 1 - (betahat1_star.smle / betahat1)
      # save to dataframe
      sim_result$betahat1_star.smle[iter] = betahat1_star.smle
      sim_result$R.model.smle[iter] = R.model.smle
      
      ########################################
      ### Output to csv ######
      ########################################
      write.csv(sim_result, 
                paste0("output/vary_prop/sims_vary_N1-",n1,"_N0-",n0,"_seed-",sim_seed,".csv"), 
                row.names = FALSE)
    }
    
    
  }
  
  toc()
}








