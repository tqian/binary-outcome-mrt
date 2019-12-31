####################################
# Created by Tianchen Qian, 2018/8/7
#
# Estimators for binary paper.
#
# Estimators include:
# - EIF
# - WCLS
# - log linear GEE
# - brm (by Richardson et al. 2017 JASA)
# Also implemented small sample correction.
#
# In this code, the ordering of parameters is (alpha, beta),
# where len(alpha) = q, and len(beta) = p.


####################################
# Update by Tianchen Qian, 2019/12/19
#
# 1. Fixed a bug in passing the result of small sample adjusted standard error estimator.


####################################
# Update by Tianchen Qian, 2019/11/28
#
# 1. Added efficient_ee_twostep, a two-step implementation for efficient_ee.
#    The first step is finding initial values of the parameters by assuming both the treatment effect model and the nuisance model is correct.
#    The second step is updating the treatment effect parameter (where the optimal weight K_t is fixed).


####################################
# Update by Tianchen Qian, 2019/3/29
#
# 1. Added GEE with exchangeable correlation structure

####################################
# Update by Tianchen Qian, 2018/9/25
#
# 1. Revised the code to allow treatment indicator = NA, when availability = 0

####################################
# Update by Tianchen Qian, 2018/8/24
#
# 1. Updated the code to incorporate availability indicator.
#
# Note of code logic in handling availability:
# the availability indicator is always multiplied with the residual r := Y - exp(alpha*Z + beta*X)
# in estimating function and in asymptotic variance calculation







library(rootSolve) # for solver function multiroot()
library(geepack) # for fitting GEE using package
library(brm) # for Richardson et al. 2017 JASA method brm()

get_alpha_beta_from_multiroot_result <- function(root, p, q)
{
    if (p == 1) {
        beta_root <- root$root[q+1]
    } else {
        beta_root <- as.matrix(root$root[(q+1) : (q+p)])
    }
    if (q == 1) {
        alpha_root <- root$root[1]
    } else {
        alpha_root <- as.matrix(root$root[1:q])
    }
    return(list(alpha = alpha_root, beta = beta_root))
}

find_change_location <- function(v){
    n <- length(v)
    if (n <= 1) {
        stop("The vector need to have length > 1.")
    }
    return(c(1, 1 + which(v[1:(n-1)] != v[2:n])))
}
# examples
# v <- c("a", "a", "b", "c", "c"); find_change_location(v)
# [1] 1 3 4


efficient_ee <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    avail_varname = NULL,
    estimator_initial_value = NULL
)
{
    ### 1. preparation ###

    sample_size <- length(unique(dta[, id_varname]))
    total_person_decisionpoint <- nrow(dta)
    
    A <- dta[, treatment_varname]
    p_t <- dta[, rand_prob_varname]
    cA <- A - p_t # centered A
    Y <- dta[, outcome_varname]
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
    
    if (is.null(avail_varname)) {
        avail <- rep(1, total_person_decisionpoint)
    } else {
        avail <- dta[, avail_varname]
    }
    
    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha
    
    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)
    
    ### 2. estimation ###
    
    estimating_equation <- function(theta) {
        alpha <- as.matrix(theta[1:q])
        beta <- as.matrix(theta[(q+1):(q+p)])
        
        exp_Zdm_alpha <- exp(Zdm %*% alpha)
        exp_Xdm_beta <- exp(Xdm %*% beta)
        exp_AXdm_beta <- exp(A * (Xdm %*% beta))
        exp_negAXdm_beta <- exp_AXdm_beta^(-1)
        exp_negXdm_beta <- exp_Xdm_beta^(-1)
        
        residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
        weight <- exp_negAXdm_beta / ( (1 - exp_Zdm_alpha) * p_t + (exp_negXdm_beta - exp_Zdm_alpha) * (1 - p_t) )
        
        ef <- rep(NA, length(theta)) # value of estimating function
        for (i in 1:q) {
            ef[i] <- sum( weight * residual * avail * Zdm[, i])
        }
        for (i in 1:p) {
            ef[q + i] <- sum( weight * residual * avail * cA * Xdm[, i])
        }
        
        ef <- ef / sample_size
        return(ef)
    }
    
    if (is.null(estimator_initial_value)) {
        estimator_initial_value <- rep(0, length = p + q)
    }

    # browser()
    solution <- tryCatch(
        {
            multiroot(estimating_equation, estimator_initial_value)
        },
        error = function(cond) {
            message("\nCatched error in multiroot inside efficient_ee():")
            message(cond)
            return(list(root = rep(NaN, p + q), msg = cond,
                        f.root = rep(NaN, p + q)))
        })
    
    estimator <- get_alpha_beta_from_multiroot_result(solution, p, q)
    alpha_hat <- as.vector(estimator$alpha)
    beta_hat <- as.vector(estimator$beta)
    
    ### 3. asymptotic variance ###
    
    ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###
    
    Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p+q, p+q))
    # Mn_summand is \frac{\partial D^{(t),T}}{\partial \theta^T} r^(t) + D^{(t),T} \frac{\partial r^(t)}{\partial \theta^T}
    # See note 2018.08.06 about small sample correction
    
    r_term_collected <- rep(NA, total_person_decisionpoint)
    D_term_collected <- matrix(NA, nrow = p+q, ncol = total_person_decisionpoint)
    partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p+q)
    
    for (it in 1:total_person_decisionpoint) {
        # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
        if (p == 1) {
            Xbeta <- Xdm[it, ] * beta_hat
        } else {
            Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
        }
        if (q == 1) {
            Zalpha <- Zdm[it, ] * alpha_hat
        } else {
            Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
        }
        
        #  old version, seems incorrect
        # W1 <- ( (1 - exp(Zalpha)) * p_t[it] + (exp(-Xbeta) - exp(Zalpha)) * (1 - p_t[it]) ) ^ (-2)
        # W2 <- exp(Zalpha - A[it] * Xbeta)
        # W3 <- - ( exp(A[it] * Xbeta) * p_t[it] * A[it] - exp(Zalpha + A[it] * Xbeta) * A[it] 
        #           + exp((A[it] - 1) * Xbeta) * (1 - p_t[it]) * (A[it] - 1) )
        
        denom <- (1 - exp(Zalpha)) * p_t[it] + (exp(-Xbeta) - exp(Zalpha)) * (1 - p_t[it])
        W1 <- exp(- A[it] * Xbeta) / (denom^2)
        W2 <- exp(Zalpha)
        W3 <- exp(-Xbeta) * (1 - p_t[it]) - A[it] * denom
        
        
        # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
        partialD_partialtheta <- matrix(NA, nrow = p + q, ncol = p + q)
        partialD_partialtheta[1:q, 1:q] <- W1 * W2 * (Zdm[it, ] %o% Zdm[it, ])
        partialD_partialtheta[1:q, (q+1):(q+p)] <- W1 * W3 * (Zdm[it, ] %o% Xdm[it, ])
        partialD_partialtheta[(q+1):(q+p), 1:q] <- W1 * W2 * cA[it] * (Xdm[it, ] %o% Zdm[it, ])
        partialD_partialtheta[(q+1):(q+p), (q+1):(q+p)] <- W1 * W3 * cA[it] * (Xdm[it, ] %o% Xdm[it, ])
        
        # r_term = r^(t) (scalar)
        r_term <- (Y[it] - exp(Zalpha + A[it] * Xbeta)) * avail[it]
        r_term_collected[it] <- r_term
        
        # D_term = D^{(t),T} (vector of length (p+q))
        D_term <- exp(- A[it] * Xbeta) / ( (1 - exp(Zalpha)) * p_t[it] + (exp(-Xbeta) - exp(Zalpha)) * (1 - p_t[it]) ) *
            c(Zdm[it, ], cA[it] * Xdm[it, ])
        D_term_collected[, it] <- D_term
        
        # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T} (vector of length (p+q))
        partialr_partialtheta <- - exp(Zalpha + A[it] * Xbeta) * c(Zdm[it, ], A[it] * Xdm[it, ]) * avail[it]
        partialr_partialtheta_collected[it, ] <- partialr_partialtheta
        
        Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
    }
    Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
    Mn_inv <- solve(Mn)
    
    ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###
    
    Sigman_summand <- array(NA, dim = c(sample_size, p+q, p+q))
    # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
    # See note 2018.08.06 about small sample correction
    
    person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)
    
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        
        Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
    }
    Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
    
    varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
    alpha_se <- sqrt(diag(varcov)[1:q])
    beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])
    
    
    ### 4. small sample correction ###
    
    Sigman_tilde <- 0
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        partialr_partialtheta_i <- partialr_partialtheta_collected[person_first_index[i] : (person_first_index[i+1] - 1), ]
        H_ii <- partialr_partialtheta_i %*% Mn_inv %*% D_term_i / sample_size
        Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii)
        
        Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
    }
    Sigman_tilde <- Sigman_tilde / sample_size
    
    varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
    alpha_se_adjusted <- sqrt(diag(varcov_adjusted)[1:q])
    beta_se_adjusted <- sqrt(diag(varcov_adjusted)[(q+1):(q+p)])
    
    ### 6. return the result with variable names ###
    
    names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames
    
    return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
                beta_se = beta_se, alpha_se = alpha_se,
                beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
                # test_result_t = test_result_t,
                # test_result_f = test_result_f,
                varcov = varcov,
                varcov_adjusted = varcov_adjusted,
                dims = list(p = p, q = q),
                f.root = solution$f.root))
}



efficient_ee_twostep <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    avail_varname = NULL
)
{
    ### 1. preparation ###
    
    sample_size <- length(unique(dta[, id_varname]))
    total_person_decisionpoint <- nrow(dta)
    
    A <- dta[, treatment_varname]
    p_t <- dta[, rand_prob_varname]
    cA <- A - p_t # centered A
    Y <- dta[, outcome_varname]
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
    
    if (is.null(avail_varname)) {
        avail <- rep(1, total_person_decisionpoint)
    } else {
        avail <- dta[, avail_varname]
    }
    
    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha
    
    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)
    
    ### 1.5. get initial values of beta and alpha to be used in the efficient estimating equation ###
    
    # This is by assuming E(Y_{t+1} | H_t, A_t) = exp(Zdm %*% alpha + A_t * Xdm %*% beta)
    # So the parameters can be found using glm(family = binomial(link = "log"))
    
    # browser()
    
    # glm_formula <- as.formula(paste0(outcome_varname, "~", paste(control_varname, collapse = "+"), "+", treatment_varname, "*(",
    #       paste(moderator_varname, collapse = "+"), ")"))
    # glm_fit <- glm(glm_formula, family = binomial(link = "log"), data = dta)
    
    # glm_fit <- glm.fit(x = cbind(Zdm, A * Xdm), y = Y, family = binomial(link = "log"))
    # geese.fit(x = cbind(Zdm, A * Xdm), y = Y, id = dta[, id_varname], family = binomial(link = "log"))
    
    # alphabeta_init <- glm_fit$coefficients # initial values of alpha and beta
    # alpha_init <- alphabeta_init[1:q]
    # beta_init <- alphabeta_init[(q+1):(q+p)]
    
    
    ## First, use an estimating equation without weights (from the derivative) to obtain the initial parameters

    ee_init <- function(theta) {
        alpha <- as.matrix(theta[1:q])
        beta <- as.matrix(theta[(q+1):(q+p)])
        
        exp_Zdm_alpha <- exp(Zdm %*% alpha)
        exp_AXdm_beta <- exp(A * Xdm %*% beta)
        
        # first time starting value through the EE without weight (from the derivative)
        residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
            
        ef <- rep(NA, length(theta)) # value of estimating function
        for (i in 1:q) {
            ef[i] <- sum( residual * avail * Zdm[, i])
        }
        for (i in 1:p) {
            ef[q + i] <- sum( residual * avail * A * Xdm[, i])
        }
        
        ef <- ef / sample_size
        return(ef)
    }
    
    solution_init <- tryCatch(
        {
            multiroot(ee_init, rep(0, p + q)) # initial value is all zero's
        },
        error = function(cond) {
            message("\nCatched error in multiroot inside efficient_ee():")
            message(cond)
            return(list(root = rep(NaN, p+q), msg = cond,
                        f.root = rep(NaN, p+q)))
        })
    
    estimator_init <- get_alpha_beta_from_multiroot_result(solution_init, p, q)
    alpha_init <- as.vector(estimator_init$alpha)
    beta_init <- as.vector(estimator_init$beta)
    
    
    ## Second, use the score equation from MLE

    ee_mle <- function(theta) {
        alpha <- as.matrix(theta[1:q])
        beta <- as.matrix(theta[(q+1):(q+p)])
        
        exp_Zdm_alpha <- exp(Zdm %*% alpha)
        exp_AXdm_beta <- exp(A * Xdm %*% beta)
        
        residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
        weight <- 1 / (1 - exp_Zdm_alpha * exp_AXdm_beta)
        
        ef <- rep(NA, length(theta)) # value of estimating function
        for (i in 1:q) {
            ef[i] <- sum( weight * residual * avail * Zdm[, i])
        }
        for (i in 1:p) {
            ef[q + i] <- sum( weight * residual * avail * A * Xdm[, i])
        }
        
        ef <- ef / sample_size
        return(ef)
    }
    
    solution_mle <- tryCatch(
        {
            multiroot(ee_mle, c(alpha_init, beta_init)) # initial value is from ee_init
        },
        error = function(cond) {
            message("\nCatched error in multiroot inside efficient_ee():")
            message(cond)
            return(list(root = rep(NaN, p+q), msg = cond,
                        f.root = rep(NaN, p+q),
                        iter = NaN, estim.precis = NaN))
        })
    
    # summary(1 / (1 - exp(Zdm %*% alpha_init) * exp(A * Xdm %*% beta_init)))
    
    if (!is.nan(solution_mle$estim.precis) & solution_mle$estim.precis < 1e-2) {
        estimator_mle <- get_alpha_beta_from_multiroot_result(solution_mle, p, q)
        alpha_mle <- as.vector(estimator_mle$alpha)
        beta_mle <- as.vector(estimator_mle$beta)
    } else {
        alpha_mle <- alpha_init
        beta_mle <- beta_init
    }

    ### 2. estimation ###
    
    exp_Zdm_alpha <- exp(Zdm %*% alpha_mle)
    exp_Xdm_beta <- exp(Xdm %*% beta_mle)
    exp_negXdm_beta <- exp_Xdm_beta^(-1)
    weight <- ( (1 - exp_Zdm_alpha) * p_t + (exp_negXdm_beta - exp_Zdm_alpha) * (1 - p_t) )^(-1)
    
    estimating_equation <- function(beta) {
        
        # only the blipping-down part uses the new beta; all other parts uses old beta
        residual <- exp( - A * (Xdm %*% beta) ) * Y - exp_Zdm_alpha
        
        ef <- rep(NA, length(beta)) # value of estimating function
        for (i in 1:p) {
            ef[i] <- sum( weight * residual * avail * cA * Xdm[, i])
        }
        
        ef <- ef / sample_size
        return(ef)
    }
    
    solution <- tryCatch(
        {
            multiroot(estimating_equation, beta_mle) # initial value here is from ee_mle
        },
        error = function(cond) {
            message("\nCatched error in multiroot inside efficient_ee():")
            message(cond)
            return(list(root = rep(NaN, p), msg = cond,
                        f.root = rep(NaN, p),
                        iter = NaN, estim.precis = NaN))
        })
    
    if (!is.nan(solution$estim.precis) & solution$estim.precis < 1e-2) {
        alpha_hat <- alpha_mle
        beta_hat <- solution$root
    } else {
        alpha_hat <- alpha_mle
        beta_hat <- beta_mle
    }

    
    ### 3. asymptotic variance ###
    
    ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###
    
    Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p+q, p+q))
    # Mn_summand is \frac{\partial D^{(t),T}}{\partial \theta^T} r^(t) + D^{(t),T} \frac{\partial r^(t)}{\partial \theta^T}
    # See note 2018.08.06 about small sample correction
    
    r_term_collected <- rep(NA, total_person_decisionpoint)
    D_term_collected <- matrix(NA, nrow = p+q, ncol = total_person_decisionpoint)
    partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p+q)
    
    for (it in 1:total_person_decisionpoint) {
        # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
        if (p == 1) {
            Xbeta <- Xdm[it, ] * beta_hat
        } else {
            Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
        }
        if (q == 1) {
            Zalpha <- Zdm[it, ] * alpha_hat
        } else {
            Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
        }
        
        #  old version, seems incorrect
        # W1 <- ( (1 - exp(Zalpha)) * p_t[it] + (exp(-Xbeta) - exp(Zalpha)) * (1 - p_t[it]) ) ^ (-2)
        # W2 <- exp(Zalpha - A[it] * Xbeta)
        # W3 <- - ( exp(A[it] * Xbeta) * p_t[it] * A[it] - exp(Zalpha + A[it] * Xbeta) * A[it] 
        #           + exp((A[it] - 1) * Xbeta) * (1 - p_t[it]) * (A[it] - 1) )
        
        denom <- (1 - exp(Zalpha)) * p_t[it] + (exp(-Xbeta) - exp(Zalpha)) * (1 - p_t[it])
        W1 <- exp(- A[it] * Xbeta) / (denom^2)
        W2 <- exp(Zalpha)
        W3 <- exp(-Xbeta) * (1 - p_t[it]) - A[it] * denom
        
        
        # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
        partialD_partialtheta <- matrix(NA, nrow = p + q, ncol = p + q)
        partialD_partialtheta[1:q, 1:q] <- W1 * W2 * (Zdm[it, ] %o% Zdm[it, ])
        partialD_partialtheta[1:q, (q+1):(q+p)] <- W1 * W3 * (Zdm[it, ] %o% Xdm[it, ])
        partialD_partialtheta[(q+1):(q+p), 1:q] <- W1 * W2 * cA[it] * (Xdm[it, ] %o% Zdm[it, ])
        partialD_partialtheta[(q+1):(q+p), (q+1):(q+p)] <- W1 * W3 * cA[it] * (Xdm[it, ] %o% Xdm[it, ])
        
        # r_term = r^(t) (scalar)
        r_term <- (Y[it] - exp(Zalpha + A[it] * Xbeta)) * avail[it]
        r_term_collected[it] <- r_term
        
        # D_term = D^{(t),T} (vector of length (p+q))
        D_term <- exp(- A[it] * Xbeta) / ( (1 - exp(Zalpha)) * p_t[it] + (exp(-Xbeta) - exp(Zalpha)) * (1 - p_t[it]) ) *
            c(Zdm[it, ], cA[it] * Xdm[it, ])
        D_term_collected[, it] <- D_term
        
        # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T} (vector of length (p+q))
        partialr_partialtheta <- - exp(Zalpha + A[it] * Xbeta) * c(Zdm[it, ], A[it] * Xdm[it, ]) * avail[it]
        partialr_partialtheta_collected[it, ] <- partialr_partialtheta
        
        Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
    }
    Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
    Mn_inv <- solve(Mn)
    
    ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###
    
    Sigman_summand <- array(NA, dim = c(sample_size, p+q, p+q))
    # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
    # See note 2018.08.06 about small sample correction
    
    person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)
    
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        
        Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
    }
    Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
    
    varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
    alpha_se <- sqrt(diag(varcov)[1:q])
    beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])
    
    
    ### 4. small sample correction ###
    
    Sigman_tilde <- 0
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        partialr_partialtheta_i <- partialr_partialtheta_collected[person_first_index[i] : (person_first_index[i+1] - 1), ]
        H_ii <- partialr_partialtheta_i %*% Mn_inv %*% D_term_i / sample_size
        Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii)
        
        Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
    }
    Sigman_tilde <- Sigman_tilde / sample_size
    
    varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
    alpha_se_adjusted <- sqrt(diag(varcov_adjusted)[1:q])
    beta_se_adjusted <- sqrt(diag(varcov_adjusted)[(q+1):(q+p)])
    
    ### 6. return the result with variable names ###
    
    names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames
    
    return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
                beta_se = beta_se, alpha_se = alpha_se,
                beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
                # test_result_t = test_result_t,
                # test_result_f = test_result_f,
                varcov = varcov,
                varcov_adjusted = varcov_adjusted,
                dims = list(p = p, q = q),
                f.root = solution$f.root))
}





efficient_ee_modified_weight <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    avail_varname = NULL,
    estimator_initial_value = NULL,
    weight_threshold = 0.95
)
{
    # weight in EE is modified to make it more stable
    
    lambda1 <- weight_threshold
    lambda2 <- weight_threshold
    
    ### 1. preparation ###
    
    sample_size <- length(unique(dta[, id_varname]))
    total_person_decisionpoint <- nrow(dta)
    
    A <- dta[, treatment_varname]
    p_t <- dta[, rand_prob_varname]
    cA <- A - p_t # centered A
    Y <- dta[, outcome_varname]
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
    
    if (is.null(avail_varname)) {
        avail <- rep(1, total_person_decisionpoint)
    } else {
        avail <- dta[, avail_varname]
    }
    
    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha
    
    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)
    
    ### 2. estimation ###
    
    estimating_equation <- function(theta) {
        alpha <- as.matrix(theta[1:q])
        beta <- as.matrix(theta[(q+1):(q+p)])
        
        exp_Zdm_alpha <- exp(Zdm %*% alpha)
        exp_Xdm_beta <- exp(Xdm %*% beta)
        exp_AXdm_beta <- exp(A * (Xdm %*% beta))
        exp_negAXdm_beta <- exp_AXdm_beta^(-1)
        exp_negXdm_beta <- exp_Xdm_beta^(-1)
        
        residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
        
        # modifying weights
        
        xi1_vector <- ifelse(exp_Zdm_alpha <= lambda1, 1 - exp_Zdm_alpha, 1 - lambda1)
        xi2_vector <- ifelse(exp_Zdm_alpha * exp_Xdm_beta <= lambda2, 1 - exp_Zdm_alpha * exp_Xdm_beta, 1 - lambda2)
        
        weight <- exp_negAXdm_beta / ( xi1_vector * p_t + exp_negXdm_beta * xi2_vector * (1 - p_t) )
        
        ef <- rep(NA, length(theta)) # value of estimating function
        for (i in 1:q) {
            ef[i] <- sum( weight * residual * avail * Zdm[, i])
        }
        for (i in 1:p) {
            ef[q + i] <- sum( weight * residual * avail * cA * Xdm[, i])
        }
        
        ef <- ef / sample_size
        return(ef)
    }
    
    if (is.null(estimator_initial_value)) {
        estimator_initial_value <- rep(0, length = p + q)
    }
    
    solution <- tryCatch(
        {
            multiroot(estimating_equation, estimator_initial_value)
        },
        error = function(cond) {
            message("\nCatched error in multiroot inside efficient_ee():")
            message(cond)
            return(list(root = rep(NaN, p + q), msg = cond,
                        f.root = rep(NaN, p + q)))
        })
    
    estimator <- get_alpha_beta_from_multiroot_result(solution, p, q)
    alpha_hat <- as.vector(estimator$alpha)
    beta_hat <- as.vector(estimator$beta)
    
    ### 3. asymptotic variance ###
    
    ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###
    
    Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p+q, p+q))
    # Mn_summand is \frac{\partial D^{(t),T}}{\partial \theta^T} r^(t) + D^{(t),T} \frac{\partial r^(t)}{\partial \theta^T}
    # See note 2018.08.06 about small sample correction
    
    r_term_collected <- rep(NA, total_person_decisionpoint)
    D_term_collected <- matrix(NA, nrow = p+q, ncol = total_person_decisionpoint)
    partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p+q)
    
    for (it in 1:total_person_decisionpoint) {
        # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
        if (p == 1) {
            Xbeta <- Xdm[it, ] * beta_hat
        } else {
            Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
        }
        if (q == 1) {
            Zalpha <- Zdm[it, ] * alpha_hat
        } else {
            Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
        }
        
        eta1 <- as.numeric(exp(Zalpha) <= lambda1)
        xi1 <- ifelse(eta1, 1 - exp(Zalpha), 1 - lambda1)
        eta2 <- as.numeric(exp(Zalpha + Xbeta) <= lambda2)
        xi2 <- ifelse(eta2, 1 - exp(Zalpha + Xbeta), 1 - lambda2)
        
        denom <- xi1 * p_t[it] + exp(-Xbeta) * xi2 * (1 - p_t[it])
        
        W1 <- exp(- A[it] * Xbeta) / (denom^2)
        W2 <- p_t[it] * eta1 + (1 - p_t[it]) * eta2
        W3 <- (1 - p_t[it]) * (exp(-Xbeta) * xi2 + exp(Zalpha) * eta2) - A[it] * denom
        
        # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
        partialD_partialtheta <- matrix(NA, nrow = p + q, ncol = p + q)
        partialD_partialtheta[1:q, 1:q] <- W1 * W2 * (Zdm[it, ] %o% Zdm[it, ])
        partialD_partialtheta[1:q, (q+1):(q+p)] <- W1 * W3 * (Zdm[it, ] %o% Xdm[it, ])
        partialD_partialtheta[(q+1):(q+p), 1:q] <- W1 * W2 * cA[it] * (Xdm[it, ] %o% Zdm[it, ])
        partialD_partialtheta[(q+1):(q+p), (q+1):(q+p)] <- W1 * W3 * cA[it] * (Xdm[it, ] %o% Xdm[it, ])
        
        # r_term = r^(t) (scalar)
        r_term <- (Y[it] - exp(Zalpha + A[it] * Xbeta)) * avail[it]
        r_term_collected[it] <- r_term
        
        # D_term = D^{(t),T} (vector of length (p+q))
        D_term <- exp(- A[it] * Xbeta) / denom * c(Zdm[it, ], cA[it] * Xdm[it, ])
        D_term_collected[, it] <- D_term
        
        # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T} (vector of length (p+q))
        partialr_partialtheta <- - exp(Zalpha + A[it] * Xbeta) * c(Zdm[it, ], A[it] * Xdm[it, ]) * avail[it]
        partialr_partialtheta_collected[it, ] <- partialr_partialtheta
        
        Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
    }
    Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
    Mn_inv <- solve(Mn)
    
    ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###
    
    Sigman_summand <- array(NA, dim = c(sample_size, p+q, p+q))
    # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
    # See note 2018.08.06 about small sample correction
    
    person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)
    
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        
        Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
    }
    Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
    
    varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
    alpha_se <- sqrt(diag(varcov)[1:q])
    beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])
    
    
    ### 4. small sample correction ###
    
    Sigman_tilde <- 0
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        partialr_partialtheta_i <- partialr_partialtheta_collected[person_first_index[i] : (person_first_index[i+1] - 1), ]
        H_ii <- partialr_partialtheta_i %*% Mn_inv %*% D_term_i / sample_size
        Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii)
        
        Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
    }
    Sigman_tilde <- Sigman_tilde / sample_size
    
    varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
    alpha_se_adjusted <- sqrt(diag(varcov_adjusted)[1:q])
    beta_se_adjusted <- sqrt(diag(varcov_adjusted)[(q+1):(q+p)])
    
    ### 6. return the result with variable names ###
    
    names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames
    
    return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
                beta_se = beta_se, alpha_se = alpha_se,
                beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
                # test_result_t = test_result_t,
                # test_result_f = test_result_f,
                varcov = varcov,
                varcov_adjusted = varcov_adjusted,
                dims = list(p = p, q = q),
                f.root = solution$f.root))
}




weighted_centered_least_square <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    avail_varname = NULL,
    rand_prob_tilde_varname = NULL, # \tilde{p}_t(1|H_t) in WCLS (variable name in the data set)
    rand_prob_tilde = NULL,         # \tilde{p}_t(1|H_t) in WCLS (numeric number or vector)
    estimator_initial_value = NULL
)
{
    ############## description ###############
    ##
    ## This function estimates the moderated treatment effect for binary outcome,
    ## and provides variance estimate.
    ##
    ## It incorporates two methods for small sample correction:
    ## 1) the usage of "Hat" matrix in the variance estimate (as in Mancl & DeRouen 2001)
    ## 2) the usage of t-distribution or F-distribution critical value with corrected degrees of freedom
    ##    (as in Liao et al. 2015)
    
    ############## arguments ###############
    ##
    ## dta...................the data set in long format
    ## id_varname............variable name for subject id (to distinguish between subjects in dta)
    ## decision_time_varname.....variable name for decision points in study
    ## outcome_varname.......variable name for outcome variable
    ## control_varname.......vector of variable names used to reduce noise (Z in the model),
    ##                       could be NULL (no control covariates)
    ## moderator_varname.....vector of variable names as effect modifiers (X in the model),
    ##                       could be NULL (no effect modifier)
    ## rand_prob_varname.....variable name for randomization probability (a column in dta)
    ## avail_varname.........variable name for availability (a column in dta)
    ##                       default to NULL (available at all decision points)
    ## rand_prob_tilde_varname.....variable name for \tilde{p}_t(1|H_t) (a column in dta)
    ##                             this is the arbitrary weight used in WCLS
    ##                             default to NULL (in which case \tilde{p}_t(1|H_t) is set to 0.5)
    ## rand_prob_tilde.............a numeric vector of the same length as dta
    ##                             this is another way to specify \tilde{p}_t(1|H_t)
    ##                             default to NULL (in which case \tilde{p}_t(1|H_t) is set to 0.5)
    ## estimator_initial_value.....initial value for the estimator in the root finding algorithm
    ##                             length is len(control_varname) + len(moderator_varname) + 2
    ##                             default to NULL (in which case the intial value = all 0's)
    
    ############## return value ###############
    ##
    ## This function returns a list of the following components:
    ##
    ## beta_hat..............estimated beta
    ## alpha_hat.............estimated alpha
    ## beta_se...............standard error for beta_hat
    ## alpha_se..............standard error for alpha_hat
    ## beta_se_adjusted......standard error for beta_hat, with small sample correction (hat matrix)
    ## alpha_se_adjusted.....standard error for alpha_hat, with small sample correction (hat matrix)
    ## varcov................estimated variance-covariance matrix for (alpha_hat, beta_hat)
    ## varcov_adjusted.......estimated variance-covariance matrix for (alpha_hat, beta_hat), with small sample correction (hat matrix)
    ## f.root................value of the estimating function at (alpha_hat, beta_hat)
    
    ##########################################
    
    ### 1. preparation ###
    
    sample_size <- length(unique(dta[, id_varname]))
    total_person_decisionpoint <- nrow(dta)
    
    if (is.null(avail_varname)) {
        avail <- rep(1, total_person_decisionpoint)
    } else {
        avail <- dta[, avail_varname]
    }
    
    A <- dta[, treatment_varname]
    # checking for NA in treatment indicator    
    if (any(is.na(A[avail == 1]))) {
        stop("Treatment indicator is NA where availability = 1.")
    }
    A[avail == 0] <- 0
    
    p_t <- dta[, rand_prob_varname]
    cA <- A - p_t # centered A
    Y <- dta[, outcome_varname]
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
    
    if (is.null(rand_prob_tilde_varname) & is.null(rand_prob_tilde)) {
        p_t_tilde <- rep(0.5, nrow(dta))
    } else if (is.null(rand_prob_tilde_varname)) {
        if (length(rand_prob_tilde) == 1) {
            p_t_tilde <- rep(rand_prob_tilde, total_person_decisionpoint)
        } else if (length(rand_prob_tilde) == total_person_decisionpoint) {
            p_t_tilde <- rand_prob_tilde
        } else {
            stop("rand_prob_tilde is of incorrect length.")
        }
    } else {
        p_t_tilde <- dta[, rand_prob_tilde_varname]
    }
    cA_tilde <- A - p_t_tilde
    
    WCLS_weight <- ifelse(A, p_t_tilde / p_t, (1 - p_t_tilde) / (1 - p_t))
    
    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha
    
    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)
    
    ### 2. estimation ###
    
    estimating_equation <- function(theta) {
        alpha <- as.matrix(theta[1:q])
        beta <- as.matrix(theta[(q+1):(q+p)])
        
        exp_Zdm_alpha <- exp(Zdm %*% alpha)
        exp_AXdm_beta <- exp(A * (Xdm %*% beta))
        
        residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
        weight <- exp_AXdm_beta^(-1)

        ef <- rep(NA, length(theta)) # value of estimating function
        for (i in 1:q) {
            ef[i] <- sum( weight * residual * avail * WCLS_weight * Zdm[, i])
        }
        for (i in 1:p) {
            ef[q + i] <- sum( weight * residual * avail * WCLS_weight * cA_tilde * Xdm[, i])
        }
        
        ef <- ef / sample_size
        return(ef)
    }
    
    if (is.null(estimator_initial_value)) {
        estimator_initial_value <- rep(0, length = p + q)
    }
    
    solution <- tryCatch(
        {
            multiroot(estimating_equation, estimator_initial_value)
        },
        error = function(cond) {
            message("\nCatched error in multiroot inside weighted_centered_least_square():")
            message(cond)
            return(list(root = rep(NaN, p + q), msg = cond,
                        f.root = rep(NaN, p + q)))
        })
    
    estimator <- get_alpha_beta_from_multiroot_result(solution, p, q)
    alpha_hat <- as.vector(estimator$alpha)
    beta_hat <- as.vector(estimator$beta)
    
    ### 3. asymptotic variance ###
    
    ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###
    
    Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p+q, p+q))
    # Mn_summand is \frac{\partial D^{(t),T}}{\partial \theta^T} r^(t) + D^{(t),T} \frac{\partial r^(t)}{\partial \theta^T}
    # See note 2018.08.06 about small sample correction
    
    r_term_collected <- rep(NA, total_person_decisionpoint)
    D_term_collected <- matrix(NA, nrow = p+q, ncol = total_person_decisionpoint)
    partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p+q)
    
    for (it in 1:total_person_decisionpoint) {
        # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
        if (p == 1) {
            Xbeta <- Xdm[it, ] * beta_hat
        } else {
            Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
        }
        if (q == 1) {
            Zalpha <- Zdm[it, ] * alpha_hat
        } else {
            Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
        }
        
        pre_multiplier <- exp(- A[it] * Xbeta) * WCLS_weight[it]
        
        # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
        partialD_partialtheta <- matrix(NA, nrow = p + q, ncol = p + q)
        partialD_partialtheta[1:q, 1:q] <- 0
        partialD_partialtheta[1:q, (q+1):(q+p)] <- - pre_multiplier * A[it] * (Zdm[it, ] %o% Xdm[it, ])
        partialD_partialtheta[(q+1):(q+p), 1:q] <- 0
        partialD_partialtheta[(q+1):(q+p), (q+1):(q+p)] <- - pre_multiplier * A[it] * cA_tilde[it] * (Xdm[it, ] %o% Xdm[it, ])
        
        # r_term = r^(t) (scalar)
        r_term <- (Y[it] - exp(Zalpha + A[it] * Xbeta)) * avail[it]
        r_term_collected[it] <- r_term
        
        # D_term = D^{(t),T} (dim = (p+q) * 1)
        D_term <- pre_multiplier * c(Zdm[it, ], cA_tilde[it] * Xdm[it, ])
        D_term_collected[, it] <- D_term
        
        # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
        partialr_partialtheta <- - exp(Zalpha + A[it] * Xbeta) * c(Zdm[it, ], A[it] * Xdm[it, ]) * avail[it]
        partialr_partialtheta_collected[it, ] <- partialr_partialtheta
        
        Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
    }
    Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
    Mn_inv <- solve(Mn)
    
    ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###
    
    Sigman_summand <- array(NA, dim = c(sample_size, p+q, p+q))
    # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
    # See note 2018.08.06 about small sample correction
    
    person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)
    
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        
        Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
    }
    Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
    
    varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
    alpha_se <- sqrt(diag(varcov)[1:q])
    beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])
    
    
    ### 4. small sample correction ###
    
    Sigman_tilde <- 0
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        partialr_partialtheta_i <- partialr_partialtheta_collected[person_first_index[i] : (person_first_index[i+1] - 1), ]
        H_ii <- partialr_partialtheta_i %*% Mn_inv %*% D_term_i / sample_size
        Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii)
        
        Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
    }
    Sigman_tilde <- Sigman_tilde / sample_size
    
    varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
    alpha_se_adjusted <- sqrt(diag(varcov_adjusted)[1:q])
    beta_se_adjusted <- sqrt(diag(varcov_adjusted)[(q+1):(q+p)])
    
    
    ### 6. return the result with variable names ###
    
    names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames
    
    return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
                beta_se = beta_se, alpha_se = alpha_se,
                beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
                varcov = varcov,
                varcov_adjusted = varcov_adjusted,
                dims = list(p = p, q = q),
                f.root = solution$f.root))
}




weighted_centered_least_square2 <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    avail_varname = NULL,
    rand_prob_tilde_varname = NULL, # \tilde{p}_t(1|H_t) in WCLS (variable name in the data set)
    rand_prob_tilde = NULL,         # \tilde{p}_t(1|H_t) in WCLS (numeric number or vector)
    estimator_initial_value = NULL
)
{
    ### 1. preparation ###
    
    sample_size <- length(unique(dta[, id_varname]))
    total_person_decisionpoint <- nrow(dta)
    
    A <- dta[, treatment_varname]
    p_t <- dta[, rand_prob_varname]
    cA <- A - p_t # centered A
    Y <- dta[, outcome_varname]
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
    
    if (is.null(avail_varname)) {
        avail <- rep(1, total_person_decisionpoint)
    } else {
        avail <- dta[, avail_varname]
    }
    
    if (is.null(rand_prob_tilde_varname) & is.null(rand_prob_tilde)) {
        p_t_tilde <- rep(0.5, nrow(dta))
    } else if (is.null(rand_prob_tilde_varname)) {
        if (length(rand_prob_tilde) == 1) {
            p_t_tilde <- rep(rand_prob_tilde, total_person_decisionpoint)
        } else if (length(rand_prob_tilde) == total_person_decisionpoint) {
            p_t_tilde <- rand_prob_tilde
        } else {
            stop("rand_prob_tilde is of incorrect length.")
        }
    } else {
        p_t_tilde <- dta[, rand_prob_tilde_varname]
    }
    cA_tilde <- A - p_t_tilde
    
    WCLS_weight <- ifelse(A, p_t_tilde / p_t, (1 - p_t_tilde) / (1 - p_t))
    
    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha
    
    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)
    
    ### 2. estimation ###
    
    estimating_equation <- function(theta) {
        alpha <- as.matrix(theta[1:q])
        beta <- as.matrix(theta[(q+1):(q+p)])
        
        exp_Zdm_alpha <- exp(Zdm %*% alpha)
        exp_AXdm_beta <- exp(A * (Xdm %*% beta))
        
        residual <- Y - exp_Zdm_alpha * exp_AXdm_beta
        weight <- exp_AXdm_beta^(-1)
        
        ef <- rep(NA, length(theta)) # value of estimating function
        for (i in 1:q) {
            ef[i] <- sum( weight * residual * avail * WCLS_weight * exp_Zdm_alpha * Zdm[, i])
        }
        for (i in 1:p) {
            ef[q + i] <- sum( weight * residual * avail * WCLS_weight * cA_tilde * Xdm[, i])
        }
        
        ef <- ef / sample_size
        return(ef)
    }
    
    if (is.null(estimator_initial_value)) {
        estimator_initial_value <- rep(0, length = p + q)
    }

    solution <- tryCatch(
        {
            multiroot(estimating_equation, estimator_initial_value)
        },
        error = function(cond) {
            message("\nCatched error in multiroot inside weighted_centered_least_square2():")
            message(cond)
            return(list(root = rep(NaN, p + q), msg = cond,
                        f.root = rep(NaN, p + q)))
        })
    
    estimator <- get_alpha_beta_from_multiroot_result(solution, p, q)
    alpha_hat <- as.vector(estimator$alpha)
    beta_hat <- as.vector(estimator$beta)
    
    ### 3. asymptotic variance ###
    
    ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###
    
    Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p+q, p+q))
    # Mn_summand is \frac{\partial D^{(t),T}}{\partial \theta^T} r^(t) + D^{(t),T} \frac{\partial r^(t)}{\partial \theta^T}
    # See note 2018.08.06 about small sample correction
    
    r_term_collected <- rep(NA, total_person_decisionpoint)
    D_term_collected <- matrix(NA, nrow = p+q, ncol = total_person_decisionpoint)
    partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p+q)
    
    for (it in 1:total_person_decisionpoint) {
        # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
        if (p == 1) {
            Xbeta <- Xdm[it, ] * beta_hat
        } else {
            Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
        }
        if (q == 1) {
            Zalpha <- Zdm[it, ] * alpha_hat
        } else {
            Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
        }
        
        pre_multiplier <- exp(- A[it] * Xbeta) * WCLS_weight[it]
        
        # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
        partialD_partialtheta <- matrix(NA, nrow = p + q, ncol = p + q)
        partialD_partialtheta[1:q, 1:q] <- pre_multiplier * exp(Zalpha) * (Zdm[it, ] %o% Zdm[it, ])
        partialD_partialtheta[1:q, (q+1):(q+p)] <- - pre_multiplier * exp(Zalpha) * A[it] * (Zdm[it, ] %o% Xdm[it, ])
        partialD_partialtheta[(q+1):(q+p), 1:q] <- 0
        partialD_partialtheta[(q+1):(q+p), (q+1):(q+p)] <- - pre_multiplier * A[it] * cA_tilde[it] * (Xdm[it, ] %o% Xdm[it, ])
        
        # r_term = r^(t) (scalar)
        r_term <- (Y[it] - exp(Zalpha + A[it] * Xbeta)) * avail[it]
        r_term_collected[it] <- r_term
        
        # D_term = D^{(t),T} (dim = (p+q) * 1)
        D_term <- pre_multiplier * c(exp(Zalpha) * Zdm[it, ], cA_tilde[it] * Xdm[it, ])
        D_term_collected[, it] <- D_term
        
        # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
        partialr_partialtheta <- - exp(Zalpha + A[it] * Xbeta) * c(Zdm[it, ], A[it] * Xdm[it, ]) * avail[it]
        partialr_partialtheta_collected[it, ] <- partialr_partialtheta
        
        Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
    }
    Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
    Mn_inv <- solve(Mn)
    
    ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###
    
    Sigman_summand <- array(NA, dim = c(sample_size, p+q, p+q))
    # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
    # See note 2018.08.06 about small sample correction
    
    person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)
    
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        
        Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
    }
    Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
    
    varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
    alpha_se <- sqrt(diag(varcov)[1:q])
    beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])
    
    
    ### 4. small sample correction ###
    
    Sigman_tilde <- 0
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        partialr_partialtheta_i <- partialr_partialtheta_collected[person_first_index[i] : (person_first_index[i+1] - 1), ]
        H_ii <- partialr_partialtheta_i %*% Mn_inv %*% D_term_i / sample_size
        Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii)
        
        Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
    }
    Sigman_tilde <- Sigman_tilde / sample_size
    
    varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
    alpha_se_adjusted <- sqrt(diag(varcov_adjusted)[1:q])
    beta_se_adjusted <- sqrt(diag(varcov_adjusted)[(q+1):(q+p)])
    
    
    ### 6. return the result with variable names ###
    
    names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames
    
    return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
                beta_se = beta_se, alpha_se = alpha_se,
                beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
                # test_result_t = test_result_t,
                # test_result_f = test_result_f,
                varcov = varcov,
                varcov_adjusted = varcov_adjusted,
                dims = list(p = p, q = q),
                f.root = solution$f.root))
}






log_linear_GEE <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    avail_varname = NULL,
    estimator_initial_value = NULL
)
{
    ### 1. preparation ###
    
    sample_size <- length(unique(dta[, id_varname]))
    total_person_decisionpoint <- nrow(dta)
    
    A <- dta[, treatment_varname]
    Y <- dta[, outcome_varname]
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
    
    if (is.null(avail_varname)) {
        avail <- rep(1, total_person_decisionpoint)
    } else {
        avail <- dta[, avail_varname]
    }
    
    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha
    
    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)
    
    ### 2. estimation ###
    
    estimating_equation <- function(theta) {
        alpha <- as.matrix(theta[1:q])
        beta <- as.matrix(theta[(q+1):(q+p)])
        
        exp_Zdm_alpha <- exp(Zdm %*% alpha)
        exp_AXdm_beta <- exp(A * (Xdm %*% beta))
        residual <- Y - exp_AXdm_beta * exp_Zdm_alpha
        weight <- (1 - exp_AXdm_beta * exp_Zdm_alpha)^(-1)
        
        ef <- rep(NA, length(theta)) # value of estimating function
        for (i in 1:q) {
            ef[i] <- sum( weight * residual * avail * Zdm[, i])
        }
        for (i in 1:p) {
            ef[q + i] <- sum( weight * residual * avail * A * Xdm[, i])
        }
        
        ef <- ef / sample_size
        return(ef)
    }
    
    if (is.null(estimator_initial_value)) {
        estimator_initial_value <- rep(0, length = p + q)
    }
    
    solution <- tryCatch(
        {
            multiroot(estimating_equation, estimator_initial_value)
        },
        error = function(cond) {
            message("\nCatched error in multiroot inside log_linear_GEE():")
            message(cond)
            return(list(root = rep(NaN, p + q), msg = cond,
                        f.root = rep(NaN, p + q)))
        })
    
    estimator <- get_alpha_beta_from_multiroot_result(solution, p, q)
    alpha_hat <- as.vector(estimator$alpha)
    beta_hat <- as.vector(estimator$beta)
    
    
    ### 3. asymptotic variance ###
    
    ### 3.1 Compute M_n matrix (M_n is the empirical expectation of the derivative of the estimating function) ###
    
    Mn_summand <- array(NA, dim = c(total_person_decisionpoint, p+q, p+q))
    # Mn_summand is \frac{\partial D^{(t),T}}{\partial \theta^T} r^(t) + D^{(t),T} \frac{\partial r^(t)}{\partial \theta^T}
    # See note 2018.08.06 about small sample correction
    
    r_term_collected <- rep(NA, total_person_decisionpoint)
    D_term_collected <- matrix(NA, nrow = p+q, ncol = total_person_decisionpoint)
    partialr_partialtheta_collected <- matrix(NA, nrow = total_person_decisionpoint, ncol = p+q)
    
    for (it in 1:total_person_decisionpoint) {
        # this is to make R code consistent whether X_it, Z_it contains more entries or is just the intercept.        
        if (p == 1) {
            Xbeta <- Xdm[it, ] * beta_hat
        } else {
            Xbeta <- as.numeric(Xdm[it, ] %*% beta_hat)
        }
        if (q == 1) {
            Zalpha <- Zdm[it, ] * alpha_hat
        } else {
            Zalpha <- as.numeric(Zdm[it, ] %*% alpha_hat)
        }
        
        pre_multiplier <- (1 - exp(Zalpha + A[it] * Xbeta))^(-2) * exp(Zalpha + A[it] * Xbeta)
        
        # partialD_partialtheta = \frac{\partial D^{(t),T}}{\partial \theta^T}, matrix of dim (p+q)*(p+q)
        partialD_partialtheta <- matrix(NA, nrow = p + q, ncol = p + q)
        partialD_partialtheta[1:q, 1:q] <- pre_multiplier * (Zdm[it, ] %o% Zdm[it, ])
        partialD_partialtheta[1:q, (q+1):(q+p)] <- pre_multiplier * A[it] * (Zdm[it, ] %o% Xdm[it, ])
        partialD_partialtheta[(q+1):(q+p), 1:q] <- pre_multiplier * A[it] * (Xdm[it, ] %o% Zdm[it, ])
        partialD_partialtheta[(q+1):(q+p), (q+1):(q+p)] <- pre_multiplier * A[it] * (Xdm[it, ] %o% Xdm[it, ])
        
        # r_term = r^(t) (scalar)
        r_term <- (Y[it] - exp(Zalpha + A[it] * Xbeta)) * avail[it]
        r_term_collected[it] <- r_term
        
        # D_term = D^{(t),T} (dim = (p+q) * 1)
        D_term <- (1 - exp(Zalpha + A[it] * Xbeta))^(-1) * c(Zdm[it, ], A[it] * Xdm[it, ])
        D_term_collected[, it] <- D_term
        
        # partialr_partialtheta = \frac{\partial r^(t)}{\partial \theta^T}
        partialr_partialtheta <- - exp(Zalpha + A[it] * Xbeta) * c(Zdm[it, ], A[it] * Xdm[it, ]) * avail[it]
        partialr_partialtheta_collected[it, ] <- partialr_partialtheta
        
        Mn_summand[it, , ] <- partialD_partialtheta * r_term + D_term %o% partialr_partialtheta
    }
    Mn <- apply(Mn_summand, c(2,3), sum) / sample_size
    Mn_inv <- solve(Mn)
    
    ### 3.2 Compute \Sigma_n matrix (\Sigma_n is the empirical variance of the estimating function) ###
    
    Sigman_summand <- array(NA, dim = c(sample_size, p+q, p+q))
    # Sigman_summand is  \sum_{t=1}^T ( D^{(t),T} r^(t) )^{\otimes 2}
    # See note 2018.08.06 about small sample correction
    
    person_first_index <- c(find_change_location(dta[, id_varname]), total_person_decisionpoint + 1)
    
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        
        Sigman_summand[i, , ] <- D_term_i %*% r_term_i %*% t(r_term_i) %*% t(D_term_i)
    }
    Sigman <- apply(Sigman_summand, c(2,3), sum) / sample_size
    
    varcov <- Mn_inv %*% Sigman %*% t(Mn_inv) / sample_size
    alpha_se <- sqrt(diag(varcov)[1:q])
    beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])
    
    
    ### 4. small sample correction ###
    
    Sigman_tilde <- 0
    for (i in 1:sample_size) {
        D_term_i <- D_term_collected[, person_first_index[i] : (person_first_index[i+1] - 1)]
        r_term_i <- matrix(r_term_collected[person_first_index[i] : (person_first_index[i+1] - 1)], ncol = 1)
        partialr_partialtheta_i <- partialr_partialtheta_collected[person_first_index[i] : (person_first_index[i+1] - 1), ]
        H_ii <- partialr_partialtheta_i %*% Mn_inv %*% D_term_i / sample_size
        Ii_minus_Hii_inv <- solve(diag(nrow(H_ii)) - H_ii)
        
        Sigman_tilde <- Sigman_tilde + D_term_i %*% Ii_minus_Hii_inv %*% r_term_i %*% t(r_term_i) %*% t(Ii_minus_Hii_inv) %*% t(D_term_i)
    }
    Sigman_tilde <- Sigman_tilde / sample_size
    
    varcov_adjusted <- Mn_inv %*% Sigman_tilde %*% t(Mn_inv) / sample_size
    alpha_se_adjusted <- sqrt(diag(varcov_adjusted)[1:q])
    beta_se_adjusted <- sqrt(diag(varcov_adjusted)[(q+1):(q+p)])
    
    
    ### 6. return the result with variable names ###
    
    names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames
    
    return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
                beta_se = beta_se, alpha_se = alpha_se,
                beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
                # test_result_t = test_result_t,
                # test_result_f = test_result_f,
                varcov = varcov,
                varcov_adjusted = varcov_adjusted,
                dims = list(p = p, q = q),
                f.root = solution$f.root))
}



log_linear_GEE_geepack <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    estimator_initial_value = NULL,
    corstr = "independence" # could also use, e.g., "exchangeable"
){
    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha
    
    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)
    
    control_summed <- paste0(control_varname, collapse = " + ")
    if (p > 1) {
        moderator_summed <- paste0("* (", paste0(moderator_varname, collapse = " + "), ")")
    } else {
        moderator_summed <- ""
    }
    
    gee_formula <- as.formula(paste(outcome_varname, "~", control_summed, "+", treatment_varname, moderator_summed))
    fit_geepack <- relRisk(gee_formula, data = dta, corstr = corstr, id = dta[, id_varname])
    
    alpha_hat <- fit_geepack$beta[1:q]
    beta_hat <- fit_geepack$beta[(q+1):(q+p)]
    varcov <- fit_geepack$vbeta
    alpha_se <- sqrt(diag(varcov)[1:q])
    beta_se <- sqrt(diag(varcov)[(q+1):(q+p)])
    alpha_se_adjusted <- alpha_se
    beta_se_adjusted <- beta_se
    
    names(alpha_hat) <- names(alpha_se) <- names(alpha_se_adjusted) <- Znames
    names(beta_hat) <- names(beta_se) <- names(beta_se_adjusted) <- Xnames
    
    return(list(beta_hat = beta_hat, alpha_hat = alpha_hat,
                beta_se = beta_se, alpha_se = alpha_se,
                beta_se_adjusted = beta_se_adjusted, alpha_se_adjusted = alpha_se_adjusted,
                # test_result_t = test_result_t,
                # test_result_f = test_result_f,
                varcov = varcov,
                dims = list(p = p, q = q),
                f.root = rep(-1, p+q)))
}




brm_MLE <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    estimator_initial_value = NULL
) {
    
    sample_size <- length(unique(dta[, id_varname]))
    total_person_decisionpoint <- nrow(dta)
    
    A <- dta[, treatment_varname]
    Y <- dta[, outcome_varname]
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
    
    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha
    
    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)
    
    est <- brm(Y, A, va = Xdm, vb = Zdm, param = "RR")
    
    return(list(beta_hat = est$point.est[1:p], alpha_hat = est$point.est[(p+1):(p+q)],
                beta_se = est$se.est[1:p], alpha_se = est$se.est[(p+1):(p+q)],
                beta_se_adjusted = est$se.est[1:p], alpha_se_adjusted = est$se.est[(p+1):(p+q)],
                # test_result_t = test_result_t,
                # test_result_f = test_result_f,
                varcov = NA,
                # varcov_ssa = asymp_varcov_ssa / sample_size,
                dims = list(p = p, q = q),
                f.root = rep(-1, p+q)))
}

brm_DR <- function(
    dta,
    id_varname,
    decision_time_varname,
    treatment_varname,
    outcome_varname,
    control_varname,
    moderator_varname,
    rand_prob_varname,
    estimator_initial_value = NULL
) {
    
    sample_size <- length(unique(dta[, id_varname]))
    total_person_decisionpoint <- nrow(dta)
    
    A <- dta[, treatment_varname]
    Y <- dta[, outcome_varname]
    Xdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, moderator_varname] ) ) # X (moderator) design matrix, intercept added
    Zdm <- as.matrix( cbind( rep(1, nrow(dta)), dta[, control_varname] ) ) # Z (control) design matrix, intercept added
    
    p <- length(moderator_varname) + 1 # dimension of beta
    q <- length(control_varname) + 1 # dimension of alpha
    
    Xnames <- c("Intercept", moderator_varname)
    Znames <- c("Intercept", control_varname)
    
    est <- brm(Y, A, va = Xdm, vb = Zdm, param = "RR", est.method = "DR")
    
    return(list(beta_hat = est$point.est[1:p], alpha_hat = est$point.est[(p+1):(p+q)],
                beta_se = est$se.est[1:p], alpha_se = est$se.est[(p+1):(p+q)],
                beta_se_adjusted = est$se.est[1:p], alpha_se_adjusted = est$se.est[(p+1):(p+q)],
                # test_result_t = test_result_t,
                # test_result_f = test_result_f,
                varcov = NA,
                # varcov_ssa = asymp_varcov_ssa / sample_size,
                dims = list(p = p, q = q),
                f.root = rep(-1, p+q)))
}
