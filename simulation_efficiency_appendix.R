rm(list = ls())

expit <- function(x){
    return(exp(x)/(1+exp(x)))
}

dgm_binary_uniform_St <- function(sample_size, total_T, rand_prob) {
    # same DGM as dgm_binary above, but faster
    
    baseline_Y_A0 <- 0.3
    
    beta_0 <- log(1/3)
    beta_1 <- 2 * log(3)
    
    df_names <- c("userid", "day", "Y", "A", "S", "prob_Y", "prob_Y_A0", "prob_A")
    
    dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
    names(dta) <- df_names
    
    dta$userid <- rep(1:sample_size, each = total_T)
    dta$day <- rep(1:total_T, times = sample_size)
    
    for (t in 1:total_T) {
        # row index for the rows corresponding to day t for every subject
        row_index <- seq(from = t, by = total_T, length = sample_size)
        
        
        dta$S[row_index] <- runif(sample_size)
        dta$prob_A[row_index] <- rep(rand_prob, sample_size)
        dta$A[row_index] <- rbinom(sample_size, 1, dta$prob_A[row_index])
        dta$prob_Y_A0[row_index] <- 0.3
        dta$prob_Y[row_index] <- dta$prob_Y_A0[row_index] * exp(dta$A[row_index] * (beta_0 + beta_1 * dta$S[row_index]))
        dta$Y[row_index] <- rbinom(sample_size, 1, dta$prob_Y[row_index])
    }
    
    return(dta)
}

beta_true <- c(log(1/3), 2 * log(3))

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

data_generating_process <- dgm_binary_uniform_St

control_vars <- "S"
moderator_vars <- "S"

nsim <- 20

rand_probs <- seq(from = 0.3, to = 0.7, by = 0.1)
sample_sizes <- c(30, 100, 200)
total_Ts <- c(10, 30, 50)

design <- expand.grid(sample_sizes, total_Ts, rand_probs)
design <- rbind(design, 
                expand.grid())
names(design) <- c("sample_size", "total_T", "rand_prob")

result_df_beta0_collected <- result_df_beta1_collected <- data.frame()

for (i_design in 1:nrow(design)) {
    print(i_design)
    
    sample_size <- design$sample_size[i_design]
    total_T <- design$total_T[i_design]
    rand_prob <- design$rand_prob[i_design]

    set.seed(123)
    
    writeLines(c(""), "~/Downloads/log.txt")
    sink("~/Downloads/log.txt", append=FALSE)
    result <- foreach(isim = 1:nsim, .combine = "c") %dorng% {
        if (isim %% 10 == 0) {
            cat(paste("Starting iteration",isim,"\n"))
        }
        dta <- data_generating_process(sample_size, total_T, rand_prob = rand_prob)
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
    
    result <- compute_result_beta(beta_true, beta, beta_se, beta_se_adjusted, moderator_vars, control_vars, significance_level = 0.05)
    result_df_beta0 <- data.frame(ss = rep(sample_size, num_estimator),
                                  total_T = rep(total_T, num_estimator),
                                  rand_prob = rep(rand_prob, num_estimator),
                                  est = ee_names,
                                  bias = result$bias[, "Intercept"],
                                  sd = result$sd[, "Intercept"],
                                  rmse = result$rmse[, "Intercept"],
                                  cp.unadj = result$coverage_prob[, "Intercept"],
                                  cp.adj = result$coverage_prob_adjusted[, "Intercept"])
    result_df_beta1 <- data.frame(ss = rep(sample_size, num_estimator),
                                  total_T = rep(total_T, num_estimator),
                                  rand_prob = rep(rand_prob, num_estimator),
                                  est = ee_names,
                                  bias = result$bias[, "S"],
                                  sd = result$sd[, "S"],
                                  rmse = result$rmse[, "S"],
                                  cp.unadj = result$coverage_prob[, "S"],
                                  cp.adj = result$coverage_prob_adjusted[, "S"])
    names(result_df_beta0) <- names(result_df_beta1) <- c("ss", "total_T", "rand_prob", "est", "bias", "sd", "rmse", "cp.unadj", "cp.adj")
    rownames(result_df_beta0) <- rownames(result_df_beta1) <- NULL
    
    result_df_beta0_collected <- rbind(result_df_beta0_collected, result_df_beta0)
    result_df_beta1_collected <- rbind(result_df_beta1_collected, result_df_beta1)

}

rm(list = ls())

load("simulation_20181227_evaluate efficiency gain along one dim submodel.rda")

design$RE_beta0 <- design$RE_beta1 <- NA
for (i_design in 1:nrow(design)) {
    design$RE_beta0[i_design] <- (result_df_beta0_collected$sd[2 * i_design - 1] / result_df_beta0_collected$sd[2 * i_design])^2
    design$RE_beta1[i_design] <- (result_df_beta1_collected$sd[2 * i_design - 1] / result_df_beta1_collected$sd[2 * i_design])^2
}

design <- design[-1,]

design_1 <- design

load("simulation_20181227_evaluate efficiency gain along one dim submodel-additional.rda")

design$RE_beta0 <- design$RE_beta1 <- NA
for (i_design in 1:nrow(design)) {
    design$RE_beta0[i_design] <- (result_df_beta0_collected$sd[2 * i_design - 1] / result_df_beta0_collected$sd[2 * i_design])^2
    design$RE_beta1[i_design] <- (result_df_beta1_collected$sd[2 * i_design - 1] / result_df_beta1_collected$sd[2 * i_design])^2
}

design <- design[-1, ]
design_2 <- design

design <- rbind(design_1, design_2)

library(ggplot2)
library(reshape2)

# p <- ggplot(design, aes(x = rand_prob, y = RE_beta1)) + geom_point()
# p + facet_grid(vars(sample_size), vars(total_T))

design_ggplot <- melt(design, id.vars = c("sample_size", "total_T", "rand_prob"))

total_T_labeller <- function(string){
    return(paste0("T = ", string))
}
sample_size_labeller <- function(string){
    return(paste0("n = ", string))
}



p <- ggplot(design_ggplot, aes(x = rand_prob, y = value, color = variable)) + geom_line()
p + facet_grid(sample_size ~ total_T,
               labeller = labeller(sample_size = sample_size_labeller,
                                   total_T = total_T_labeller)) +
    scale_x_continuous(breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
    scale_color_manual(name  = "estimand",
                       labels = c(expression(paste(psi[1])),
                                  expression(paste(psi[0]))),
                       values = c("#000000", "#E69F00")) +
    xlab("randomization probability") +
    ylab("relative efficiency") +
    theme_bw()

ggsave("simulation_20181227_evaluate efficiency gain along one dim submodel.png",
       width = 6, height = 4)
