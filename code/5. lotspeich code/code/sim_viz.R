###########################################
### MCAR
###########################################
setwd("code/5. lotspeich code")
library(latex2exp)
library(ggplot2)

###########################################
# MCAR 
###########################################

#############################
## Weight model with intercept only

sim_res = read.csv("mcar_weight1_sim_res.csv")

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
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param"), 
                                labels = c("GS (NP)", "GS (P)", 
                                           "CC (NP)", "CC (P)", 
                                           "IPW (NP)", "IPW (P)",
                                           "SMLE (P)")), 
                parametric = !grepl(pattern = "(NP)", 
                                    x = method))
jpeg('mcar_weight1.jpg')
# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white")) + 
  ggtitle(label = "Boxplot of estimates under MCAR with intercept weight")
dev.off()

#############################
## Weight model with Y only

sim_res = read.csv("mcar_weightY_sim_res.csv")

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
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param"), 
                                labels = c("GS (NP)", "GS (P)", 
                                           "CC (NP)", "CC (P)", 
                                           "IPW (NP)", "IPW (P)",
                                           "SMLE (P)")), 
                parametric = !grepl(pattern = "(NP)", 
                                    x = method))
jpeg('mcar_weightY.jpg')
# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white")) + 
  ggtitle(label = "Boxplot of estimates under MCAR with Y weight model")
dev.off()

#############################
## Weight model with Z only

sim_res = read.csv("mcar_weightZ_sim_res.csv")

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
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param"), 
                                labels = c("GS (NP)", "GS (P)", 
                                           "CC (NP)", "CC (P)", 
                                           "IPW (NP)", "IPW (P)",
                                           "SMLE (P)")), 
                parametric = !grepl(pattern = "(NP)", 
                                    x = method))
jpeg('mcar_weightZ.jpg')
# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white")) + 
  ggtitle(label = "Boxplot of estimates under MCAR with Z weight model")
dev.off()

#############################
## Weight model with Y + Z

sim_res = read.csv("mcar_weightY+Z_sim_res.csv")

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
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param"), 
                                labels = c("GS (NP)", "GS (P)", 
                                           "CC (NP)", "CC (P)", 
                                           "IPW (NP)", "IPW (P)",
                                           "SMLE (P)")), 
                parametric = !grepl(pattern = "(NP)", 
                                    x = method))
jpeg('mcar_weightY+Z.jpg')
# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white")) + 
  ggtitle(label = "Boxplot of estimates under MCAR with Y and Z weight model")
dev.off()

###########################################
### MAR given Y
###########################################

#############################
## Weight model with intercept only

sim_res = read.csv("mar_contY_weight1_sim_res.csv")

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
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param"), 
                                labels = c("GS (NP)", "GS (P)", 
                                           "CC (NP)", "CC (P)", 
                                           "IPW (NP)", "IPW (P)",
                                           "SMLE (P)")), 
                parametric = !grepl(pattern = "(NP)", 
                                    x = method))
jpeg('mar_contY_weight1.jpg')
# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white")) + 
  ggtitle(label = "Boxplot of estimates under MAR by Y using weight intercept model")
dev.off()

#############################
## Weight model with Y only

sim_res = read.csv("mar_contY_weightY_sim_res.csv")

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
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param"), 
                                labels = c("GS (NP)", "GS (P)", 
                                           "CC (NP)", "CC (P)", 
                                           "IPW (NP)", "IPW (P)",
                                           "SMLE (P)")), 
                parametric = !grepl(pattern = "(NP)", 
                                    x = method))
jpeg('mar_contY_weightY.jpg')
# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white")) + 
  ggtitle(label = "Boxplot of estimates under MAR by Y using weight Y model")
dev.off()

#############################
## Weight model with Z only

sim_res = read.csv("mar_contY_weightZ_sim_res.csv")

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
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param"), 
                                labels = c("GS (NP)", "GS (P)", 
                                           "CC (NP)", "CC (P)", 
                                           "IPW (NP)", "IPW (P)",
                                           "SMLE (P)")), 
                parametric = !grepl(pattern = "(NP)", 
                                    x = method))
jpeg('mar_contY_weightZ.jpg')
# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white")) + 
  ggtitle(label = "Boxplot of estimates under MAR by Y using weight Z model")
dev.off()

#############################
## Weight model with Y + Z

sim_res = read.csv("mar_contY_weightY+Z_sim_res.csv")

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
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param"), 
                                labels = c("GS (NP)", "GS (P)", 
                                           "CC (NP)", "CC (P)", 
                                           "IPW (NP)", "IPW (P)",
                                           "SMLE (P)")), 
                parametric = !grepl(pattern = "(NP)", 
                                    x = method))
jpeg('mar_contY_weightY+Z.jpg')
# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white")) + 
  ggtitle(label = "Boxplot of estimates under MAR by Y using weight Y and Z model")
dev.off()

###########################################
### MAR given Z
###########################################

#############################
## Weight model with intercept only

sim_res = read.csv("mar_Z_weight1_sim_res.csv")

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
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param"), 
                                labels = c("GS (NP)", "GS (P)", 
                                           "CC (NP)", "CC (P)", 
                                           "IPW (NP)", "IPW (P)",
                                           "SMLE (P)")), 
                parametric = !grepl(pattern = "(NP)", 
                                    x = method))
jpeg('mar_Z_weight1.jpg')
# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white")) + 
  ggtitle(label = "Boxplot of estimates under MAR by Z using weight intercept model")
dev.off()

#############################
## Weight model with Y only

sim_res = read.csv("mar_Z_weightY_sim_res.csv")

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
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param"), 
                                labels = c("GS (NP)", "GS (P)", 
                                           "CC (NP)", "CC (P)", 
                                           "IPW (NP)", "IPW (P)",
                                           "SMLE (P)")), 
                parametric = !grepl(pattern = "(NP)", 
                                    x = method))
jpeg('mar_Z_weightY.jpg')
# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white")) + 
  ggtitle(label = "Boxplot of estimates under MAR by Z using weight Y model")
dev.off()

#############################
## Weight model with Z only

sim_res = read.csv("mar_Z_weightZ_sim_res.csv")

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
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param"), 
                                labels = c("GS (NP)", "GS (P)", 
                                           "CC (NP)", "CC (P)", 
                                           "IPW (NP)", "IPW (P)",
                                           "SMLE (P)")), 
                parametric = !grepl(pattern = "(NP)", 
                                    x = method))
jpeg('mar_Z_weightZ.jpg')
# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white")) + 
  ggtitle(label = "Boxplot of estimates under MAR by Z using weight Z model")
dev.off()

#############################
## Weight model with Y + Z

sim_res = read.csv("mar_Z_weightY+Z_sim_res.csv")

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
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param"), 
                                labels = c("GS (NP)", "GS (P)", 
                                           "CC (NP)", "CC (P)", 
                                           "IPW (NP)", "IPW (P)",
                                           "SMLE (P)")), 
                parametric = !grepl(pattern = "(NP)", 
                                    x = method))
jpeg('mar_Z_weightY+Z.jpg')
# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white")) + 
  ggtitle(label = "Boxplot of estimates under MAR by Z using weight Y and Z model")
dev.off()
