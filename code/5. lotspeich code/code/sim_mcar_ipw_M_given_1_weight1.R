# clear environment
rm(list=ls())
# load libraries and utility codes
source("~/Research/missing_surrogate/code/5. lotspeich code/R.s.miss.R")
source("~/Research/missing_surrogate/code/5. lotspeich code/gen_data.R")

# seed
sim_seed = 30
# set num of reps
num_reps = 1000
# set sample size
n1 = 1000 ### control group size
n0 = 1000 ### treatment group size
# prop of missingness for each treatment group, from 0 to 1
prop_m0 = 0.7
prop_m1 = 0.7

sim_res = data.frame(r = 1:1000, 
                     gs_nonparam_delta = NA, gs_nonparam_delta.s = NA, gs_nonparam_R.s = NA, 
                     gs_param_delta = NA, gs_param_delta.s = NA, gs_param_R.s = NA, 
                     cc_nonparam_delta = NA, cc_nonparam_delta.s = NA, cc_nonparam_R.s = NA, 
                     cc_param_delta = NA, cc_param_delta.s = NA, cc_param_R.s = NA, 
                     ipw_nonparam_delta = NA, ipw_nonparam_delta.s = NA, ipw_nonparam_R.s = NA, 
                     ipw_param_delta = NA, ipw_param_delta.s = NA, ipw_param_R.s = NA)

for (r in 1:1000) {
  # Generate data 
  data = gen.data(n1=n1, n0=n0) 
  
  # Define vectors for outcomes/surrogates in untreated/treated 
  s1 = data$s1
  y1 = data$y1
  s0 = data$s0
  y0 = data$y0
  
  ##########################################################
  #Estimates with complete data ############################
  ##########################################################
  ## Estimate R with nonparametric approach (gold standard)
  Rnonparam = Rsurrogate::R.s.estimate(sone=s1, 
                                       szero=s0, 
                                       yone = y1, 
                                       yzero = y0, 
                                       type = "robust")
  sim_res[r, c("gs_nonparam_delta", "gs_nonparam_delta.s", "gs_nonparam_R.s")] = with(Rnonparam, 
                                                                                      c(delta, delta.s, R.s))
  
  ## Estimate R with parametric approach (gold standard) ###
  Rparam = Rsurrogate::R.s.estimate(sone=s1, 
                                    szero=s0, 
                                    yone = y1, 
                                    yzero = y0, 
                                    type = "model")
  sim_res[r, c("gs_param_delta", "gs_param_delta.s", "gs_param_R.s")] = with(Rparam, 
                                                                             c(delta, delta.s, R.s))
  
  # Simulate non-missingness indicators ###################
  ## Under MCAR, everybody has 30% missingness probability
  m1 = rbinom(n1, 1, prop_m1) 
  m0 = rbinom(n0, 1, prop_m0)
  s0[m0==0] = NA ### make them missing
  s1[m1==0] = NA ### make them missing
  
  ##########################################################
  #Estimates with incomplete data ##########################
  ##########################################################
  ## Estimate R with nonparametric approach (complete case)
  Rnonparam_miss = Rsurrogate::R.s.estimate(sone = s1[m1==1], 
                                            szero = s0[m0==1], 
                                            yone = y1[m1==1], 
                                            yzero = y0[m0==1], 
                                            type = "robust")
  sim_res[r, c("cc_nonparam_delta", "cc_nonparam_delta.s", "cc_nonparam_R.s")] = with(Rnonparam_miss, 
                                                                                      c(delta, delta.s, R.s))
  
  ## Estimate R with parametric approach (complete case)
  Rparam_miss = Rsurrogate::R.s.estimate(sone = s1[m1==1], 
                                         szero = s0[m0==1], 
                                         yone = y1[m1==1], 
                                         yzero = y0[m0==1], 
                                         type = "model")
  sim_res[r, c("cc_param_delta", "cc_param_delta.s", "cc_param_R.s")] = with(Rparam_miss, 
                                                                             c(delta, delta.s, R.s))
  ipw_dat = data.frame(m = c(m1, m0), 
                       y = c(y1, y0), 
                       z = rep(x = c(1, 0), each = 1000))
  ## Calculate weights for IPW approaches
  ipw_fit = glm(formula = m ~ 1, 
                family = "binomial", 
                data = ipw_dat)
  w1 = predict(object = ipw_fit, 
               newdata = data.frame(m = m1, 
                                    y = y1, 
                                    z = 1),
               type = "response")
  w0 = predict(object = ipw_fit, 
               newdata = data.frame(m = m0, 
                                    y = y0, 
                                    z = 0),
               type = "response")
  
  ## Estimate R with nonparametric approach (IPW)
  Rnonparam_miss_ipw = R.s.miss(sone = s1, 
                                szero = s0, 
                                yone = y1, 
                                yzero = y0, 
                                type = "robust", 
                                wone = w1, 
                                wzero = w0)
  sim_res[r, c("ipw_nonparam_delta", "ipw_nonparam_delta.s", "ipw_nonparam_R.s")] = with(Rnonparam_miss_ipw, 
                                                                                         c(delta, delta.s, R.s))
  
  ## Estimate R with parametric approach (IPW)
  Rparam_miss_ipw = R.s.miss(sone = s1, 
                             szero = s0, 
                             yone = y1, 
                             yzero = y0, 
                             type = "model", 
                             wone = w1, 
                             wzero = w0)
  sim_res[r, c("ipw_param_delta", "ipw_param_delta.s", "ipw_param_R.s")] = with(Rparam_miss_ipw,
                                                                                c(delta, delta.s, R.s))
  ## Estimate R with semiparametric approach (SMLE)
  Rparam_miss_smle = R.s.miss(sone = s1, 
                              szero = s0,
                              yone = y1,
                              yzero = y0, 
                              type = "model") 
  sim_res[r, c("smle_param_delta", "smle_param_delta.s", "smle_param_R.s")] = with(Rparam_miss_smle, c(delta, delta.s, R.s))
  
  ## Save 
  sim_res |> 
    write.csv("~/Research/missing_surrogate/code/5. lotspeich code/mcar_weight1_sim_res.csv", 
              row.names = FALSE)
}

