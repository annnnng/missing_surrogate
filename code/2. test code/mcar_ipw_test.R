# clear environment
rm(list=ls())
# load libraries and utility codes
source("code/gen_data.R")
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
sim_seed = 4
# set num of reps
num_reps = 1000
# set sample size
n1 = 1000 ### control group size
n0 = 1000 ### treatment group size
# prop of missingness for each treatment group, from 0 to 1
prop_obs_s0 = 0.3
prop_obs_s1 = 0.5

sink(paste0("output/sim_vary_prop_obs", format(Sys.time(), "%d-%b-%Y %H.%M"), ".txt"))

for (prop_obs_s0 in c(0.3, 0.5, 0.7, 1.0)) {
  
  for (prop_obs_s1 in c(0.3, 0.5, 0.7, 1.0)){
    
    print(paste("System time =", Sys.time()))
    tic(paste("Sims with prop s0 =", prop_obs_s0, "- prop s1 =", prop_obs_s1))
    
    # set seed
    set.seed(sim_seed)
    
    ### Store sims data
    sim_result = data.frame(
      sim = paste(sim_seed, 1:num_reps, sep = "-"),
      n1 = n1, n0 = n0, prop_obs_s0 = prop_obs_s0, prop_obs_s1 = prop_obs_s1,
      R.gold = NA, R.complete = NA, R.ipw.ipw = NA, R.ipw.complete = NA,
      R.smle.ipw = NA, R.smle.complete = NA
    )
    
    for(iter in 1:num_reps){
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
      
      ##########################################################
      #Estimates with incomplete data ##########################
      ##########################################################
      
      
      # Simulate non-missingness indicators ###################
      ## Under MCAR, everybody has 30% missingness probability
      m1.1 = rbinom(n1,1,prop_obs_s1) 
      m1.0 = rbinom(n0,1,prop_obs_s0)
      df$Smiss[c(m1.0, m1.1) == 0] = NA
      ## Missingness weight
      w.y1 = rep(0.7, n1)
      w.y0 = rep(0.7, n0)
      df$weights = c(w.y0,w.y1)
      
      ### Estimate R with nonparametric approach (complete case)
      R_miss_complete = Rsurrogate::R.s.estimate(sone=s1[m1.1==1], 
                                                 szero=s0[m1.0==1], 
                                                 yone = y1[m1.1==1], 
                                                 yzero = y0[m1.0==1], 
                                                 type = "model")
      sim_result$R.complete[iter] = R_miss_complete$R.s
      
      ## IPW
      
      ### Estimate R with IPW only (semiparam)
      
      #### Fit Y ~ Z + S
      lm_ipw_semiparam = lm(formula = Y ~ Z*Smiss, 
                            data = df, 
                            weights = weights)
      summary(lm_ipw_semiparam)
      
      ##### Get coefficients
      gamma0 = as.numeric(lm_ipw_semiparam$coefficients[1])
      gamma1 = as.numeric(lm_ipw_semiparam$coefficients[3])
      gamma2 = gamma0 + as.numeric(lm_ipw_semiparam$coefficients[2])
      gamma3 = gamma1 + as.numeric(lm_ipw_semiparam$coefficients[4])
      
      #### Fit S ~ Z with IPW
      lm_surr_weighted = lm(formula = Smiss ~ Z, data = df, 
                            weights = weights)
      summary(lm_surr_weighted)
      ##### Get coefficients
      alpha0_weighted = as.numeric(lm_surr_weighted$coefficients[1])
      alpha1_weighted = alpha0_weighted + as.numeric(lm_surr_weighted$coefficients[2])
      
      #### Calculate PTE 
      delta = (gamma2 - gamma0) + gamma1 * (alpha1_weighted - alpha0_weighted) + (gamma3 - gamma1) * alpha1_weighted
      deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0_weighted
      sim_result$R.ipw.ipw[iter]  = 1 - deltaS / delta 
      
      #### Fit S ~ Z complete case
      lm_surr_comp = lm(formula = Smiss ~ Z, data = df)
      ##### Get coefficients
      alpha0_comp = as.numeric(coefficients(lm_surr_comp)[1])
      alpha1_comp = alpha0_comp + as.numeric(coefficients(lm_surr_comp)[2])
      
      #### Calculate PTE 
      delta = (gamma2 - gamma0) + gamma1 * (alpha1_comp - alpha0_comp) + (gamma3 - gamma1) * alpha1_comp 
      deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0_comp
      sim_result$R.ipw.complete[iter]  = 1 - deltaS / delta 
      
      ## smle
      
      ### Estimate R with smle & IPW (semiparam)
      
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
      fit$coefficients
      gamma0 = as.numeric(fit$coefficients[1])
      gamma1 = as.numeric(fit$coefficients[2])
      gamma2 = gamma0 + as.numeric(fit$coefficients[4])
      gamma3 = gamma1 + as.numeric(fit$coefficients[3])
      
      #### Calculate PTE 
      delta = (gamma2 - gamma0) + gamma1 * (alpha1_weighted - alpha0_weighted) + (gamma3 - gamma1) * alpha1_weighted
      deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0_weighted
      sim_result$R.smle.ipw[iter]  = 1 - deltaS / delta  
      
      
      #### Calculate PTE 
      delta = (gamma2 - gamma0) + gamma1 * (alpha1_comp - alpha0_comp) + (gamma3 - gamma1) * alpha1_comp 
      deltaS = (gamma2 - gamma0) + (gamma3 - gamma1) * alpha0_comp
      sim_result$R.smle.complete[iter]  = 1 - deltaS / delta 
      
    }
    
    write.csv(sim_result, 
              paste0("output/vary_prop/sims_vary_prop-s0-",prop_obs_s0,"_s1-",prop_obs_s1,"_seed-",sim_seed,".csv"), 
              row.names = FALSE)
  }
  
  toc()
}

sink()







