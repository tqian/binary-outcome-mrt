
# for binary outcome

# Tianchen Qian, 2019.03.03

# In this code, I constructed examples where:
# Z_t (covariate) is an exogenous AR(1) process
# p_t depends on Z_t through logistic regression
# E(Y_t+1 | H_t, A_t = 0) depends on Z_t, Y_t, A_t-1 through logistic regression or exp()
# treatment effect is constant

expit <- function(x){
    return(exp(x)/(1+exp(x)))
}

prob_clip <- function(vector, min, max) {
    vector[which(vector < min)] <- min
    vector[which(vector > max)] <- max
    return(vector)
}

dgm_binary_ar1_covariate <- function(sample_size, total_T, eta = 0.5, gamma = 0.4, y_control_link = c("expit", "exp")) {
    # eta determines how prob(A_t) is affected by Z_t
    # eta = 0 means A_t independent of Z_t
    # gamma determines how prob_Y_A0 depends on Z_t, gamma = 0 means no dependence
    
    y_control_link <- match.arg(y_control_link)
    
    beta_0 <- 0.1
    
    df_names <- c("userid", "day", "Y", "A", "Z", "prob_Y", "prob_Y_A0", "prob_A", "Y_lag1")
    
    dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
    names(dta) <- df_names
    
    dta$userid <- rep(1:sample_size, each = total_T)
    dta$day <- rep(1:total_T, times = sample_size)
    
    dta$Z <- as.vector(replicate(sample_size, arima.sim(model = list(ar = 0.5), n = total_T)))
    
    # first time point
    t <- 1
    row_index <- seq(from = t, by = total_T, length = sample_size)
    
    dta$prob_A[row_index] <- prob_clip(expit(eta * dta$Z[row_index]), min = 0.2, max = 0.8)
    # dta$prob_A[row_index] <- 0.5
    dta$A[row_index] <- rbinom(sample_size, 1, dta$prob_A[row_index])
    if (y_control_link == "expit") {
        dta$prob_Y_A0[row_index] <- prob_clip(expit(gamma * dta$Z[row_index]), min = 0.1, max = 0.8)
    } else if (y_control_link == "exp") {
        dta$prob_Y_A0[row_index] <- prob_clip(exp(gamma * (dta$Z[row_index] - 3)), min = 0.1, max = 0.8)
    }
    dta$prob_Y[row_index] <- dta$prob_Y_A0[row_index] * exp(dta$A[row_index] * beta_0)
    dta$Y[row_index] <- rbinom(sample_size, 1, dta$prob_Y[row_index])
    dta$Y_lag1[row_index] <- 0
    
    for (t in 2:total_T) {
        # row index for the rows corresponding to day t for every subject
        row_index <- seq(from = t, by = total_T, length = sample_size)
        row_index_pre <- seq(from = t-1, by = total_T, length = sample_size)
        
        dta$prob_A[row_index] <- prob_clip(expit(eta * dta$Z[row_index]), min = 0.2, max = 0.8)
        # dta$prob_A[row_index] <- 0.5
        dta$A[row_index] <- rbinom(sample_size, 1, dta$prob_A[row_index])
        if (y_control_link == "expit") {
            dta$prob_Y_A0[row_index] <- prob_clip(
                expit(-0.5 + gamma * dta$Z[row_index] + 0.2 * dta$Y[row_index_pre] + 0.2 * dta$A[row_index_pre]),
                min = 0.1, max = 0.8)
        } else if (y_control_link == "exp") {
            dta$prob_Y_A0[row_index] <- prob_clip(
                exp(- 0.4 + gamma * (dta$Z[row_index] - 3) + 0.2 * dta$Y[row_index_pre] + 0.2 * dta$A[row_index_pre]),
                min = 0.1, max = 0.8)
        }
        
        dta$prob_Y[row_index] <- dta$prob_Y_A0[row_index] * exp(dta$A[row_index] * beta_0)
        dta$Y[row_index] <- rbinom(sample_size, 1, dta$prob_Y[row_index])
        dta$Y_lag1[row_index] <- dta$Y[row_index_pre]
    }
    
    return(dta)
}

## true beta
beta_0_true <- 0.1

# try out the range of Y
if (0) {
    set.seed(123)
    dta <- dgm_binary_ar1_covariate(100, 30, eta = - 0.5, gamma = 0.5)
    summary(dta$prob_Y)
    summary(dta$prob_Y_A0)
    summary(dta$prob_A)
    hist(dta$prob_Y_A0, xlim = c(0,1))
    hist(dta$prob_Y, xlim = c(0,1))
    hist(dta$prob_A, xlim = c(0,1))
}

