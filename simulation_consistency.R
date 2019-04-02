# code to conduct simulation 1 and create table for binary outcome paper

# Tianchen Qian
# 2018.08.12

# simulation part is copied from simulation.R
# table creation part is copied from U-stat paper

# update on 2019.02.06: to include eif_modified_weight

# update on 2019.03.29: to include GEE with exchangeable correlation structure


##### simulation part #####

rm(list = ls())

source("estimators.R")
source("estimators_robust_adhocery.R")

compute_result_beta <- function(beta_true, beta, beta_se, beta_se_adjusted, moderator_vars, control_vars, significance_level,
                                na.rm = FALSE) {
    
    beta_true_array <- array(NA, dim = dim(beta), dimnames = dimnames(beta))
    for (ind1 in 1:dim(beta_true_array)[1]) {
        for (ind3 in 1:dim(beta_true_array)[3]) {
            beta_true_array[ind1, , ind3] <- beta_true
        }
    }
    
    p <- length(moderator_vars) + 1
    q <- length(control_vars) + 1
    
    bias <- apply(beta - beta_true_array, c(1,2), mean, na.rm = na.rm)
    sd <- apply(beta, c(1,2), sd, na.rm = na.rm)
    rmse <- apply(beta - beta_true_array, c(1,2), function(v) sqrt(mean(v^2, na.rm = na.rm)))
    
    critical_factor <- qnorm(1 - significance_level/2)
    ci_left <- beta - critical_factor * beta_se
    ci_right <- beta + critical_factor * beta_se
    coverage_prob <- apply((ci_left < beta_true_array) & (ci_right > beta_true_array),
                           c(1,2), mean, na.rm = na.rm)
    
    critical_factor_adj <- qt(1 - significance_level/2, df = sample_size - 1 - q)
    ci_left_adj <- beta - critical_factor_adj * beta_se_adjusted
    ci_right_adj <- beta + critical_factor_adj * beta_se_adjusted
    coverage_prob_adj <- apply((ci_left_adj < beta_true_array) & (ci_right_adj > beta_true_array),
                               c(1,2), mean, na.rm = na.rm)
    
    return(list(bias = bias, sd = sd, rmse = rmse, coverage_prob = coverage_prob, coverage_prob_adjusted = coverage_prob_adj))
}

source("dgm_binary_categorical_covariate.R")
data_generating_process <- dgm_binary_categorical_covariate

library(tidyverse)

library(foreach)
library(doMC)
library(doRNG)

max_cores <- 16
registerDoMC(min(detectCores() - 1, max_cores))

sample_sizes <- c(30, 50, 100)
total_T <- 30
nsim <- 1000

control_vars <- "S"
moderator_vars <- c()

result_df_collected <- data.frame()

for (i_ss in 1:length(sample_sizes)) {
    
    sample_size <- sample_sizes[i_ss]
    
    writeLines(c(""), "~/Downloads/log.txt")
    sink("~/Downloads/log.txt", append=FALSE)
    result <- foreach(isim = 1:nsim, .combine = "c") %dorng% {
        if (isim %% 10 == 0) {
            cat(paste("Starting iteration",isim,"\n"))
        }
        dta <- data_generating_process(sample_size, total_T)
        fit_wcls <- weighted_centered_least_square(
            dta = dta,
            id_varname = "userid",
            decision_time_varname = "day",
            treatment_varname = "A",
            outcome_varname = "Y",
            control_varname = control_vars,
            moderator_varname = moderator_vars,
            rand_prob_varname = "prob_A",
            rand_prob_tilde_varname = NULL,
            rand_prob_tilde = 0.2,
            estimator_initial_value = NULL
        )
        
        fit_eif <- efficient_ee(
            dta = dta,
            id_varname = "userid",
            decision_time_varname = "day",
            treatment_varname = "A",
            outcome_varname = "Y",
            control_varname = control_vars,
            moderator_varname = moderator_vars,
            rand_prob_varname = "prob_A",
            estimator_initial_value = c(fit_wcls$alpha_hat, fit_wcls$beta_hat)
        )
        
        fit_eif_modified <- efficient_ee_modified_weight(
            dta = dta,
            id_varname = "userid",
            decision_time_varname = "day",
            treatment_varname = "A",
            outcome_varname = "Y",
            control_varname = control_vars,
            moderator_varname = moderator_vars,
            rand_prob_varname = "prob_A",
            estimator_initial_value = c(fit_wcls$alpha_hat, fit_wcls$beta_hat)
        )
        
        fit_gee_ind <- log_linear_GEE_geepack(
            dta = dta,
            id_varname = "userid",
            decision_time_varname = "day",
            treatment_varname = "A",
            outcome_varname = "Y",
            control_varname = control_vars,
            moderator_varname = moderator_vars,
            estimator_initial_value = NULL,
            corstr = "independence"
        )
        
        fit_gee_exch <- log_linear_GEE_geepack(
            dta = dta,
            id_varname = "userid",
            decision_time_varname = "day",
            treatment_varname = "A",
            outcome_varname = "Y",
            control_varname = control_vars,
            moderator_varname = moderator_vars,
            estimator_initial_value = NULL,
            corstr = "exchangeable"
        )
        
        output <- list(list(fit_wcls = fit_wcls, fit_eif = fit_eif, fit_eif_modified = fit_eif_modified, fit_gee_ind = fit_gee_ind, fit_gee_exch = fit_gee_exch))
    }
    sink()
    
    ee_names <- c("wcls", "eif", "eif_modified", "gee_ind", "gee_exch")
    alpha_names <- c("Intercept", control_vars)
    beta_names <- c("Intercept", moderator_vars)
    num_estimator <- length(ee_names)
    
    
    alpha <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls$alpha_hat, l$fit_eif$alpha_hat, l$fit_eif_modified$alpha_hat, l$fit_gee_ind$alpha_hat, l$fit_gee_exch$alpha_hat),
                                                              nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, alpha_names))))
    alpha_se <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls$alpha_se, l$fit_eif$alpha_se, l$fit_eif_modified$alpha_se, l$fit_gee_ind$alpha_se, l$fit_gee_exch$alpha_se),
                                                                 nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, alpha_names))))
    alpha_se_adjusted <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls$alpha_se_adjusted, l$fit_eif$alpha_se_adjusted, l$fit_eif_modified$alpha_se_adjusted, l$fit_gee_ind$alpha_se_adjusted, l$fit_gee_exch$alpha_se_adjusted),
                                                                          nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, alpha_names))))
    beta <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls$beta_hat, l$fit_eif$beta_hat, l$fit_eif_modified$beta_hat, l$fit_gee_ind$beta_hat, l$fit_gee_exch$beta_hat),
                                                             nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, beta_names))))
    beta_se <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls$beta_se, l$fit_eif$beta_se, l$fit_eif_modified$beta_se, l$fit_gee_ind$beta_se, l$fit_gee_exch$beta_se),
                                                                nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, beta_names))))
    beta_se_adjusted <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls$beta_se_adjusted, l$fit_eif$beta_se_adjusted, l$fit_eif_modified$beta_se_adjusted, l$fit_gee_ind$beta_se_adjusted, l$fit_gee_exch$beta_se_adjusted),
                                                                         nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, beta_names))))
    
    result <- compute_result_beta(beta_true_marginal, beta, beta_se, beta_se_adjusted, moderator_vars, control_vars, significance_level = 0.05)
    result_df <- data.frame(ss = rep(sample_size, num_estimator),
                            est = ee_names,
                            bias = result$bias,
                            sd = result$sd,
                            rmse = result$rmse,
                            cp.unadj = result$coverage_prob,
                            cp.adj = result$coverage_prob_adjusted)
    names(result_df) <- c("ss", "est", "bias", "sd", "rmse", "cp.unadj", "cp.adj")
    rownames(result_df) <- NULL
    
    result_df_collected <- rbind(result_df_collected, result_df)
}

saveRDS(result_df_collected, file = "result_simulation_1.RDS")

##### create tables for paper #####

library(reshape)
library(kableExtra)
library(knitr)

result_df_collected <- readRDS("result_simulation_1.RDS")


result_df_collected <- result_df_collected[, c(2, 1, 3:ncol(result_df_collected))]
result_df_collected$est <- factor(result_df_collected$est, c("wcls", "eif", "eif_modified", "gee_ind", "gee_exch"))
result_df_collected <- result_df_collected[order(result_df_collected$est, result_df_collected$ss), ]

rownames(result_df_collected) <- NULL

result_df_collected$bias <- round(result_df_collected$bias, 3)
result_df_collected$sd <- round(result_df_collected$sd, 3)
result_df_collected$rmse <- round(result_df_collected$rmse, 3)
result_df_collected$cp.unadj <- round(result_df_collected$cp.unadj, 2)
result_df_collected$cp.adj <- round(result_df_collected$cp.adj, 2)

sink("table_generation/simulation_1.txt", append=FALSE)
mycaption <- "caption for simulation 1"
latex_code <- kable(result_df_collected, format = "latex", booktabs = T, align = "c", caption = mycaption) %>%
    # add_header_above(c("est", "sample.size", "bias", "sd", "rmse", "cp.unadj", "cp.adj")) %>%
    # column_spec(1, bold=T) %>%
    collapse_rows(columns = 1, latex_hline = "major")
print(latex_code)
sink()

# then run replace_text-boot.py to replace certain text in the output table
