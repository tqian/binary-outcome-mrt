rm(list = ls())

source("dgm_binary_ar1_covariate.R")

# try out the range of Y
if (0) {
    set.seed(123)
    dta <- dgm_binary_uniform_St(100, 30, rand_prob = 0.5)
    summary(dta$prob_Y)
    summary(dta$prob_A)
}


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





library(foreach)
library(doMC)
library(doRNG)

max_cores <- 16
registerDoMC(min(detectCores() - 1, max_cores))

source("estimators.R")

data_generating_process <- dgm_binary_ar1_covariate

# control_vars <- c("Z", "Y_lag1")
# control_vars <- c("Y_lag1")
control_vars <- c("Z")
moderator_vars <- NULL

nsim <- 1000

sample_sizes <- c(50)
# total_Ts <- c(10, 30)
total_Ts <- c(20)
# sample_sizes <- c(100)
# total_Ts <- c(30)
etas <- c(-0.5, 0, 0.5)
# etas <- -0.5
gammas <- seq(from = 0.1, to = 0.5, by = 0.1)

design <- expand.grid(sample_sizes, total_Ts, etas, gammas)
names(design) <- c("sample_size", "total_T", "eta", "gamma")

result_df_beta0_collected <- result_df_beta1_collected <- data.frame()

for (i_design in 1:nrow(design)) {
# for (i_design in 4:4) {
    print(i_design)
    
    sample_size <- design$sample_size[i_design]
    total_T <- design$total_T[i_design]
    eta <- design$eta[i_design]
    gamma <- design$gamma[i_design]

    set.seed(20190303)
    
    writeLines(c(""), "~/Downloads/log.txt")
    sink("~/Downloads/log.txt", append=FALSE)
    result <- foreach(isim = 1:nsim, .combine = "c") %dorng% {
        if (isim %% 10 == 0) {
            cat(paste("Starting iteration",isim,"\n"))
        }
        dta <- data_generating_process(sample_size, total_T, eta, gamma, "expit")
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
            rand_prob_tilde = 0.5,
            estimator_initial_value = NULL
        )
        
        # fit_eif <- efficient_ee(
        #     dta = dta,
        #     id_varname = "userid",
        #     decision_time_varname = "day",
        #     treatment_varname = "A",
        #     outcome_varname = "Y",
        #     control_varname = control_vars,
        #     moderator_varname = moderator_vars,
        #     rand_prob_varname = "prob_A",
        #     estimator_initial_value = c(fit_wcls$alpha_hat, fit_wcls$beta_hat)
        # )
        fit_eif <- fit_wcls
        
        fit_eif_modified <- efficient_ee_modified_weight(
            dta = dta,
            id_varname = "userid",
            decision_time_varname = "day",
            treatment_varname = "A",
            outcome_varname = "Y",
            control_varname = control_vars,
            moderator_varname = moderator_vars,
            rand_prob_varname = "prob_A",
            estimator_initial_value = c(fit_wcls$alpha_hat, fit_wcls$beta_hat),
            weight_threshold = 0.97
        )
        
        output <- list(list(fit_wcls = fit_wcls, fit_eif = fit_eif, fit_eif_modified = fit_eif_modified))
    }
    sink()
    
    ee_names <- c("wcls", "eif", "eif_modified")
    alpha_names <- c("Intercept", control_vars)
    beta_names <- c("Intercept", moderator_vars)
    num_estimator <- length(ee_names)
    
    alpha <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls$alpha_hat, l$fit_eif$alpha_hat, l$fit_eif_modified$alpha_hat),
                                                              nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, alpha_names))))
    alpha_se <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls$alpha_se, l$fit_eif$alpha_se, l$fit_eif_modified$alpha_se),
                                                                 nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, alpha_names))))
    alpha_se_adjusted <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls$alpha_se_adjusted, l$fit_eif$alpha_se_adjusted, l$fit_eif_modified$alpha_se_adjusted),
                                                                          nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, alpha_names))))
    beta <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls$beta_hat, l$fit_eif$beta_hat, l$fit_eif_modified$beta_hat),
                                                             nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, beta_names))))
    beta_se <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls$beta_se, l$fit_eif$beta_se, l$fit_eif_modified$beta_se),
                                                                nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, beta_names))))
    beta_se_adjusted <- simplify2array(lapply(result, function(l) matrix(c(l$fit_wcls$beta_se_adjusted, l$fit_eif$beta_se_adjusted, l$fit_modified$beta_se_adjusted),
                                                                         nrow = length(ee_names), byrow = TRUE, dimnames = list(ee_names, beta_names))))
    
    result <- compute_result_beta(beta_0_true, beta, beta_se, beta_se_adjusted, moderator_vars, control_vars, significance_level = 0.05)
    result_df_beta0 <- data.frame(ss = rep(sample_size, num_estimator),
                                  total_T = rep(total_T, num_estimator),
                                  eta = rep(eta, num_estimator),
                                  gamma = rep(gamma, num_estimator),
                                  est = ee_names,
                                  bias = result$bias[, "Intercept"],
                                  sd = result$sd[, "Intercept"],
                                  rmse = result$rmse[, "Intercept"],
                                  cp.unadj = result$coverage_prob[, "Intercept"],
                                  cp.adj = result$coverage_prob_adjusted[, "Intercept"])
    names(result_df_beta0) <- c("ss", "total_T", "eta", "gamma", "est", "bias", "sd", "rmse", "cp.unadj", "cp.adj")
    rownames(result_df_beta0) <- NULL
    
    result_df_beta0_collected <- rbind(result_df_beta0_collected, result_df_beta0)

}

result_df_beta0_collected

save.image("simulation_20190303_evaluate efficiency gain along one dim submodel(expit).rda")

rm(list = ls())

load("simulation_20190303_evaluate efficiency gain along one dim submodel(exp).rda")

design$RE_beta0 <- NA
for (i_design in 1:nrow(design)) {
    design$RE_beta0[i_design] <- (result_df_beta0_collected$sd[3 * i_design - 2] / result_df_beta0_collected$sd[3 * i_design])^2
}
design$link <- "q-exp"
design_exp <- design

load("simulation_20190303_evaluate efficiency gain along one dim submodel(expit).rda")

design$RE_beta0 <- NA
for (i_design in 1:nrow(design)) {
    design$RE_beta0[i_design] <- (result_df_beta0_collected$sd[3 * i_design - 2] / result_df_beta0_collected$sd[3 * i_design])^2
}
design$link <- "q-expit"
design_expit <- design


design <- rbind(design_exp, design_expit)

design$sample_size <- NULL
design$total_T <- NULL

design <- subset(design, gamma == 0.10 | gamma == 0.50 | (gamma > 0.29 & gamma < 0.31))

library(ggplot2)
library(reshape2)

# p <- ggplot(design, aes(x = rand_prob, y = RE_beta1)) + geom_point()
# p + facet_grid(vars(sample_size), vars(total_T))

design_ggplot <- melt(design, id.vars = c("eta", "gamma", "link"))



p <- ggplot(design, aes(x = gamma, y = RE_beta0, color = factor(eta))) + geom_line()
p + facet_grid(. ~ link) +
    scale_color_manual(name  = expression(eta),
                       values = c("#000000", "#E69F00", "#56B4E9")) +
    xlab(expression(gamma)) +
    ylab("relative efficiency") +
    theme_bw()

ggsave("simulation_20190303_evaluate efficiency gain along one dim submodel.png",
       width = 4, height = 2)

