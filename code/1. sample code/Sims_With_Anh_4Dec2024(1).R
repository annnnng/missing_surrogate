library(Rsurrogate)
data("d_example")

# Transform the data to be long 
long_d_example = data.frame(Y = c(d_example$y0, d_example$y1), 
                            Z = rep(c(0, 1), each = 500), 
                            S = c(d_example$s0.a, d_example$s1.a))

# Replicate 
## Fit Y ~ Z 
modYgivZ = lm(formula = Y ~ Z, 
              data = long_d_example)
betahat1 = coefficients(modYgivZ)[2]

## Fit Y ~ Z + S
modYgivZS = lm(formula = Y ~ Z + S, 
               data = long_d_example)
betahat1_star = coefficients(modYgivZS)[2]

## Calculate PTE = R_F = 1 - \beta_1^* / \beta_1
RF = 1 - (betahat1_star / betahat1)

# Make some of the surrogates missing 
long_d_example$Smiss = long_d_example$S ## Make a new surrogate column 
long_d_example$Smiss[1:75] = NA ## Make the first 75 / 500 (15% of the treated) missing 
long_d_example$Smiss[851:1000] = NA ## Make the last 150 / 500 (30% of the untreated) missing 

# Fit Ran's two-phase SMLE for Y ~ Z + S 
long_d_example$bs1 = as.numeric(long_d_example$Z == 0) ## Create bs1 = I(Z = 0), treated 
long_d_example$bs2 = as.numeric(long_d_example$Z == 1) ## Create bs2 = I(Z = 1), untreated
#devtools::install_github("dragontaoran/TwoPhaseReg")
fit = TwoPhaseReg::smle(Y = "Y", ## outcome  
                        X = "Smiss", ## "expensive" covariate (with missing data)
                        Z = "Z", ## "inexpensive" covariate (without missing data)
                        Bspline_Z = c("bs1", "bs2"), ## B-splines / indicators for Z 
                        data = long_d_example, ## name of the dataset 
                        noSE = FALSE, ## do you want standard error estimates? 
                        TOL = 1E-4, ## used to define convergence
                        MAX_ITER = 500, ## used to keep it from blowing up 
                        model = "linear" ## type of outcome model - we want OLS
                        )
betahat1_star_miss = fit$coefficients[3]

## Calculate PTE = R_F = 1 - \beta_1^*(miss) / \beta_1
RF_miss = 1 - (betahat1_star_miss / betahat1)
