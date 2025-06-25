setwd("code/5. lotspeich code")
library(latex2exp)
library(ggplot2)
source("boxplot.R")

###########################################
# MCAR 
###########################################

# read data for different weight model
sim_res_weight1 = read.csv("mcar_weight1_sim_res.csv")
sim_res_weight1$model = "intercept"
sim_res_weightY = read.csv("mcar_weightY_sim_res.csv")
sim_res_weightY$model = "Y"
sim_res_weightZ = read.csv("mcar_weightZ_sim_res.csv")
sim_res_weightZ$model = "Z"
sim_res_weightYZ = read.csv("mcar_weightY+Z_sim_res.csv")
sim_res_weightYZ$model = "YZ"

# combine data with a column for weight model 
sim_res = rbind(sim_res_weight1, sim_res_weightY)
sim_res = rbind(sim_res, sim_res_weightZ)
sim_res = rbind(sim_res, sim_res_weightYZ)

# parse column heading 

# Make them long 
res_long = sim_res |> 
  tidyr::pivot_longer(cols = gs_nonparam_delta:smle_param_R.s, 
                      names_to = "method_quantity", values_to = "est") |> 
  dplyr::mutate(quantity = sub(pattern = ".*_", 
                               replacement = "", 
                               x = method_quantity), 
                truth = dplyr::case_when(
                  quantity == "delta" ~ 12, 
                  quantity == "delta.s" ~ 6, 
                  .default = 0.5), 
# nice label with latex                
#                quantity = factor(x = quantity, 
#                                  levels = c("delta", "delta.s", "R.s"), 
#                                  labels = c(TeX("Quantity: $\\Delta$"), TeX("Quantity: $\\Delta_S$"), TeX("Quantity: $R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(parametric = factor(x = !grepl(pattern = "nonparam", x = method), 
                                    levels = c(FALSE, TRUE), 
                                    labels = c("PTE Estimator: Nonparametric",
                                               "PTE Estimator: Parametric")), 
                method = factor(x = method, 
                                levels = c("gs_nonparam", "cc_nonparam", "ipw_nonparam",
                                           "gs_param", "cc_param", "ipw_param", "smle_param"), 
                                labels = c("Gold Standard", "Complete Case", "IPW",
                                           "Gold Standard",  "Complete Case", "IPW",  "SMLE")),
                model = factor(x = model, 
                                levels = c("intercept", "Y", "Z", "YZ"), 
                                labels = c(TeX("Weight Model: $M \\sim 1$"), TeX("Weight Model: $M \\sim Y$"), 
                                           TeX("Weight Model: $M \\sim Z$"), TeX("Weight Model: $M \\sim Y + Z$")
                                           )
                               )
                )


jpeg('mcar_boxplot_est_model.jpg')
# make boxplot
res_long |>
  dplyr::filter(quantity == "R.s") |> 
  ggplot(aes(x = method, y = est, fill = method)) + 
  geom_hline(aes(yintercept = 0.5), linetype = 2, color = "black") + 
  geom_boxplot() + 
  xlab("Method") + 
  ylab("Estimate") + 
  facet_grid(cols = vars(parametric), 
             rows = vars(model), 
             scales = "free",
             labeller = labeller(parametric = label_value, 
                                model = label_parsed)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  scale_fill_manual(name = "Method:", values = cols, guide = "none") + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white"))
dev.off()


###########################################
# MAR on Z
###########################################

# read data for different weight model
sim_res_weight1 = read.csv("mar_Z_weight1_sim_res.csv")
sim_res_weight1$model = "intercept"
sim_res_weightY = read.csv("mar_Z_weightY_sim_res.csv")
sim_res_weightY$model = "Y"
sim_res_weightZ = read.csv("mar_Z_weightZ_sim_res.csv")
sim_res_weightZ$model = "Z"
sim_res_weightYZ = read.csv("mar_Z_weightY+Z_sim_res.csv")
sim_res_weightYZ$model = "YZ"

# combine data with a column for weight model 
sim_res = rbind(sim_res_weight1, sim_res_weightY)
sim_res = rbind(sim_res, sim_res_weightZ)
sim_res = rbind(sim_res, sim_res_weightYZ)

# parse column heading 

# Make them long 
res_long = sim_res |> 
  tidyr::pivot_longer(cols = gs_nonparam_delta:smle_param_R.s, 
                      names_to = "method_quantity", values_to = "est") |> 
  dplyr::mutate(quantity = sub(pattern = ".*_", 
                               replacement = "", 
                               x = method_quantity), 
                truth = dplyr::case_when(
                  quantity == "delta" ~ 12, 
                  quantity == "delta.s" ~ 6, 
                  .default = 0.5), 
                # nice label with latex                
                #                quantity = factor(x = quantity, 
                #                                  levels = c("delta", "delta.s", "R.s"), 
                #                                  labels = c(TeX("Quantity: $\\Delta$"), TeX("Quantity: $\\Delta_S$"), TeX("Quantity: $R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(parametric = factor(x = !grepl(pattern = "nonparam", x = method), 
                                    levels = c(FALSE, TRUE), 
                                    labels = c("PTE Estimator: Nonparametric",
                                               "PTE Estimator: Parametric")), 
                method = factor(x = method, 
                                levels = c("gs_nonparam", "cc_nonparam", "ipw_nonparam",
                                           "gs_param", "cc_param", "ipw_param", "smle_param"), 
                                labels = c("Gold Standard", "Complete Case", "IPW",
                                           "Gold Standard",  "Complete Case", "IPW",  "SMLE")),
                model = factor(x = model, 
                               levels = c("intercept", "Y", "Z", "YZ"), 
                               labels = c(TeX("Weight Model: $M \\sim 1$"), TeX("Weight Model: $M \\sim Y$"), 
                                          TeX("Weight Model: $M \\sim Z$"), TeX("Weight Model: $M \\sim Y + Z$")
                               )
                )
  )


jpeg('mar_Z_boxplot_est_model.jpg')
# make boxplot
res_long |>
  dplyr::filter(quantity == "R.s") |> 
  ggplot(aes(x = method, y = est, fill = method)) + 
  geom_hline(aes(yintercept = 0.5), linetype = 2, color = "black") + 
  geom_boxplot() + 
  xlab("Method") + 
  ylab("Estimate") + 
  facet_grid(cols = vars(parametric), 
             rows = vars(model), 
             scales = "free",
             labeller = labeller(parametric = label_value, 
                                 model = label_parsed)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  scale_fill_manual(name = "Method:", values = cols, guide = "none") + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white"))
dev.off()

###########################################
# MAR on Y
###########################################

# read data for different weight model
sim_res_weight1 = read.csv("mar_contY_weight1_sim_res.csv")
sim_res_weight1$model = "intercept"
sim_res_weightY = read.csv("mar_contY_weightY_sim_res.csv")
sim_res_weightY$model = "Y"
sim_res_weightZ = read.csv("mar_contY_weightZ_sim_res.csv")
sim_res_weightZ$model = "Z"
sim_res_weightYZ = read.csv("mar_contY_weightY+Z_sim_res.csv")
sim_res_weightYZ$model = "YZ"

# combine data with a column for weight model 
sim_res = rbind(sim_res_weight1, sim_res_weightY)
sim_res = rbind(sim_res, sim_res_weightZ)
sim_res = rbind(sim_res, sim_res_weightYZ)

# parse column heading 

# Make them long 
res_long = sim_res |> 
  tidyr::pivot_longer(cols = gs_nonparam_delta:smle_param_R.s, 
                      names_to = "method_quantity", values_to = "est") |> 
  dplyr::mutate(quantity = sub(pattern = ".*_", 
                               replacement = "", 
                               x = method_quantity), 
                truth = dplyr::case_when(
                  quantity == "delta" ~ 12, 
                  quantity == "delta.s" ~ 6, 
                  .default = 0.5), 
                # nice label with latex                
                #                quantity = factor(x = quantity, 
                #                                  levels = c("delta", "delta.s", "R.s"), 
                #                                  labels = c(TeX("Quantity: $\\Delta$"), TeX("Quantity: $\\Delta_S$"), TeX("Quantity: $R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(parametric = factor(x = !grepl(pattern = "nonparam", x = method), 
                                    levels = c(FALSE, TRUE), 
                                    labels = c("PTE Estimator: Nonparametric",
                                               "PTE Estimator: Parametric")), 
                method = factor(x = method, 
                                levels = c("gs_nonparam", "cc_nonparam", "ipw_nonparam",
                                           "gs_param", "cc_param", "ipw_param", "smle_param"), 
                                labels = c("Gold Standard", "Complete Case", "IPW",
                                           "Gold Standard",  "Complete Case", "IPW",  "SMLE")),
                model = factor(x = model, 
                               levels = c("intercept", "Y", "Z", "YZ"), 
                               labels = c(TeX("Weight Model: $M \\sim 1$"), TeX("Weight Model: $M \\sim Y$"), 
                                          TeX("Weight Model: $M \\sim Z$"), TeX("Weight Model: $M \\sim Y + Z$")
                               )
                )
  )


jpeg('mar_contY_boxplot_est_model.jpg')
# make boxplot
res_long |>
  dplyr::filter(quantity == "R.s") |> 
  ggplot(aes(x = method, y = est, fill = method)) + 
  geom_hline(aes(yintercept = 0.5), linetype = 2, color = "black") + 
  geom_boxplot() + 
  xlab("Method") + 
  ylab("Estimate") + 
  facet_grid(cols = vars(parametric), 
             rows = vars(model), 
             scales = "free",
             labeller = labeller(parametric = label_value, 
                                 model = label_parsed)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  scale_fill_manual(name = "Method:", values = cols, guide = "none") + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white"))
dev.off()


