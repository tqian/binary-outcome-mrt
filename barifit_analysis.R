# Exploratory data analysis with Barifit data
# Tianchen Qian
# 11/30/2018


# update 2020.01.19
# Revised the exploratory plots.

# update 2019.12.04
# Added two exploratory plots (tx by day, tx by user).


# update 2019.04.07
# Add code to use GEE to analyze BariFit data.
# However, cannot open the file
# data_filename <- paste0(sys.var$mbox, "Bari-Fit data/barifit_csv_files/food_track_analysis_data_correct.csv")
# this time. Error message is
#   Error in read.table(file = file, header = header, sep = sep, quote = quote,  : 
#   no lines available in input


# update 12/28/2018
# added some exploratory analysis to look at the data




##### questions to ask Pedja
# 1. What's the difference between "food_track_analysis_data_correct.csv" and "food track analysis data.csv",
#    both files are in "/Bari-Fit data/barifit_csv_files/" folder
# 2. How to get other covariates for food track analysis?


rm(list = ls())

library(tidyverse)


sys.var <- switch(Sys.info()["sysname"],
                  "Windows" = list(locale = "English",
                                   mbox = "Z:/BariFit data/"),
                  "Darwin" = list(locale = "en_US",
                                  mbox = "/Volumes/dav/BariFit data/"),
                  "Linux" = list(locale = "en_US.UTF-8",
                                 mbox = "~/mbox/BariFit data/"))


##### activity suggestion data set #####

suggest <- read.csv(paste0(sys.var$mbox, "Bari-Fit data/barifit_csv_files/MRT_activity_suggestion_data.csv"))

# Get the start date and end date of each individual.
# This is to align other time-varying covariates with the food track reminder data set (which only has day but not dates).

extract_date_suggest <- function(string) {
    date <- as.Date(substr(as.character(string), 1, 7), format = "%d%b%y")
    return(date)
}

summarize_date <- function(vector_of_dates) {
    # return: first date, last date, has_all_dates_in_between, range_including_first_and_last
    first_date <- min(vector_of_dates)
    last_date <- max(vector_of_dates)
    range_including_first_and_last <- last_date - first_date + 1
    has_all_dates_in_between <- (length(unique(vector_of_dates)) == range_including_first_and_last)
    missing_days <- range_including_first_and_last - length(unique(vector_of_dates))
    
    return(data.frame(first_date = first_date,
                      last_date = last_date,
                      range_including_first_and_last = range_including_first_and_last,
                      has_all_dates_in_between = has_all_dates_in_between,
                      missing_days = missing_days))
}

suggest$date <- as.Date(NA)
for (i in 1:nrow(suggest)) {
    suggest$date[i] <- extract_date_suggest(suggest$server_sent_dtm[i])
}

date_range <- by(suggest$date, suggest$study_id, summarize_date)
date_range <- do.call(rbind, date_range)


##### Get baseline variables #####

demog <- read.csv(paste0(sys.var$mbox, "Bari-Fit data/barifit_csv_files/barifit_demog_EMR.csv"),
                  stringsAsFactors = FALSE)
demog$DOB <- as.Date(demog$DOB, format = "%m/%d/%Y")

### covariate: age 
hist(demog$DOB, breaks = 50, freq = T) # there is no obvious seperation in age

### covariate: gener
table(demog$Gender) # 43 Female, 8 Male

### covariate: diabetes_condition
table(demog$diabetes_condition)
# BTH ICD  NP TRT 
#   8   4  30   9 

### covariate: bmi

weight_emr <- read.csv(paste0(sys.var$mbox, "Bari-Fit data/barifit_csv_files/barifit_weight_EMR.csv"),
                       stringsAsFactors = FALSE)
weight_emr$measure_date <- as.Date(weight_emr$measure_date, format = "%m/%d/%y")

user_ids <- rownames(date_range)
baseline_bmi <- data.frame(Study_ID = user_ids, bmi = NA, measure_date = as.Date(NA))
for (i in 1:length(user_ids)) {
    user_id <- user_ids[i]
    dta_subset <- subset(weight_emr, Study_ID == user_id)
    tmp <- which(dta_subset$measure_date <= date_range$first_date[i])
    baseline_measure_date_index <- tmp[length(tmp)]
    baseline_bmi$bmi[i] <- dta_subset$bmi[baseline_measure_date_index]
    baseline_bmi$measure_date[i] <- dta_subset$measure_date[baseline_measure_date_index]
}



##### prepare data set for analysis #####

data_filename <- paste0(sys.var$mbox, "Bari-Fit data/barifit_csv_files/food_track_analysis_data_correct.csv")

dta <- read.csv(data_filename, stringsAsFactors = FALSE)

# descriptive statistics

str(dta)

summary(dta)

table(dta$study_id)

### add covariate: Y lag 1

dta$Y_previous_day <- c(0, dta$Y[1:(nrow(dta) - 1)])
dta$Y_previous_day[dta$Day == 112] <- 0

### add baseline covariate: gender, age, bmi
dta$gender <- NA
dta$age <- NA
dta$bmi <- NA

user_ids <- unique(dta$study_id)

for (i in 1:length(user_ids)) {
    user_id <- user_ids[i]
    demog_row <- demog[demog$Study_ID == user_id, ]
    date_range_row <- date_range[rownames(date_range) == user_id, ]
    bmi_row <- baseline_bmi[baseline_bmi$Study_ID == user_id, ]
    row_index_in_dta <- which(dta$study_id == user_id)
    
    dta$gender[row_index_in_dta] <- demog_row$Gender
    dta$age[row_index_in_dta] <- 
        as.numeric(difftime(date_range_row$first_date, demog_row$DOB, units = "weeks")) / 52.25
    dta$bmi[row_index_in_dta] <- bmi_row$bmi
}

dta$gender <- ifelse(dta$gender == "M", 0, 1)
dta$study_id <- as.factor(dta$study_id)


dta$prob_A <- 0.5
dta$Day <- dta$Day - 1


##### check randomization #####

tmp <- by(dta, dta$study_id, function(data) mean(data$message_sent))
hist(tmp, main = "Histogram of empirical randomization \n prob by individual", xlab = "empirical rand prob")
mean(dta$message_sent)

summary(glm(message_sent ~ gender + age + bmi + Day + Y_previous_day, data = dta, family = binomial))




##### preliminary GEE analysis to select control variables #####


library(geepack)


fit <- geeglm(Y ~ Day + Y_previous_day + gender + age + bmi, data = dta,
              corstr = "independence", id = study_id, family = binomial)
summary(fit)

# Coefficients:
#                 Estimate   Std.err    Wald Pr(>|W|)    
# (Intercept)    -2.199548  1.209133   3.309  0.06889 .  
# Day            -0.010449  0.002116  24.382  7.9e-07 ***
# Y_previous_day  3.295458  0.260893 159.554  < 2e-16 ***
# gender         -1.133466  0.387793   8.543  0.00347 ** 
# age             0.019658  0.010484   3.516  0.06079 .  
# bmi             0.031457  0.025911   1.474  0.22473  



##### Appendix H: exploratory plots #####


myggfont <- function(legend_text_size = 14,
                     legend_title_size = 16,
                     axis_text_size = 16,
                     axis_title_size = 20,
                     plot_title_size = 24) {
  ff <- theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size, face="bold"),
              plot.title = element_text(size = plot_title_size))
  return(ff)
}

# average completion rate under treatment, under no treatment, and on average

mean(dta$Y[dta$message_sent == 1]) # 0.512
mean(dta$Y[dta$message_sent == 0]) # 0.5118
mean(dta$Y) # 0.5119

# plotting the average completion rate under treatment and under no treatment by Day

tmp_day <- dta %>% group_by(Day) %>% summarize(avg_Y_trt = mean(Y[message_sent == 1]), avg_Y_ctl = mean(Y[message_sent == 0]))
tmp_day$avg_Y_diff <- tmp_day$avg_Y_trt - tmp_day$avg_Y_ctl
tmp_day$tx_negative <- (tmp_day$avg_Y_diff <= 0)
tmp_day2 <- dta %>% group_by(Day, message_sent) %>% summarize(avg_Y = mean(Y))

ggplot(tmp_day) +
  geom_segment(aes(x = Day, xend = Day, y = avg_Y_trt, yend = avg_Y_ctl, linetype = tx_negative)) +
  geom_point(aes(x = Day, y = avg_Y, shape = as.factor(1 - message_sent)), data = tmp_day2, size = 2.5) +
  theme_bw() +
  ylab("Food log completion rate") +
  xlab("Day in study") +
  scale_linetype_manual(name  = "Crude treatment effect",
                        labels = c("positive", "negative"),
                        values = c(1, 2)) +
  scale_shape_manual(name  = "Average outcome under",
                     labels = c("treatment", "no treatment"),
                     values = c(1, 2)) +
  coord_cartesian(ylim = c(0,1)) +
  myggfont() +
  theme(legend.position="bottom")


ggsave("plot/avg_Y_by_day.pdf", width = 11, height = 5.5)


# plotting the average completion rate under treatment and under no treatment by user

tmp_user <- dta %>% group_by(study_id) %>% summarize(avg_Y_trt = mean(Y[message_sent == 1]), avg_Y_ctl = mean(Y[message_sent == 0]))
tmp_user$avg_Y_diff <- tmp_user$avg_Y_trt - tmp_user$avg_Y_ctl
tmp_user$rank <- rank(- tmp_user$avg_Y_diff)
tmp_user$tx_negative <- (tmp_user$avg_Y_diff <= 0)
# tmp_user2 <- dta %>% group_by(study_id, message_sent) %>% summarize(avg_Y = mean(Y))
# tmp_user2 <- tmp_user2[order(tmp_user2$avg_Y_diff, decreasing = TRUE), ]
# tmp_user2$rank <- rep(tmp_user$rank, each = 2)

ggplot(tmp_user) +
  geom_segment(aes(x = rank, xend = rank, y = avg_Y_trt, yend = avg_Y_ctl, linetype = tx_negative)) +
  geom_point(aes(x = rank, y = avg_Y, shape = as.factor(1 - message_sent)), data = tmp_user2, size = 3) +
  theme_bw() +
  ylab("Food log completion rate") +
  xlab("Participant id (sorted by crude treatment effect)") +
  scale_linetype_manual(name  = "Crude treatment effect",
                        labels = c("positive", "negative"),
                        values = c(1, 2)) +
  scale_shape_manual(name  = "Average outcome under",
                     labels = c("treatment", "no treatment"),
                     values = c(1, 2)) +
  coord_cartesian(ylim = c(0,1)) +
  myggfont() +
  theme(legend.position="bottom")

ggsave("plot/avg_Y_by_user.pdf", width = 11, height = 5.5)



##### Section 7: analysis in main paper #####

dta$prob_A <- 0.5

### analysis 1: marginal excursion effect

control_vars <- c("Day", "Y_previous_day", "gender")
moderator_vars <- NULL

fit_wcls <- weighted_centered_least_square(
  dta = dta,
  id_varname = "study_id",
  decision_time_varname = "Day",
  treatment_varname = "message_sent",
  outcome_varname = "Y",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = 0.5,
  estimator_initial_value = NULL
)

t_quantile <- qt(0.975, 45 - 1 - 2) 

# estimator, SE, 95% CI, p-value

fit_wcls$beta_hat # 0.01416758
fit_wcls$beta_se_adjusted # 0.02089641
rbind(fit_wcls$beta_hat - t_quantile * fit_wcls$beta_se_adjusted,
      fit_wcls$beta_hat + t_quantile * fit_wcls$beta_se_adjusted)
# (-0.02800308, 0.05633824)
2 * pt(abs(fit_wcls$beta_hat) / fit_wcls$beta_se_adjusted, 45 - 1 - 2, lower.tail = FALSE) # 0.5014956



### analysis 2: causal excursion effect moderation by Day


control_vars <- c("Day", "Y_previous_day", "gender")
moderator_vars <- "Day"

fit_wcls <- weighted_centered_least_square(
  dta = dta,
  id_varname = "study_id",
  decision_time_varname = "Day",
  treatment_varname = "message_sent",
  outcome_varname = "Y",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = 0.5,
  estimator_initial_value = NULL
)

t_quantile <- qt(0.975, 45 - 1 - 2) 

# estimator, SE, 95% CI, p-value
fit_wcls$beta_hat
# Intercept          Day 
# 0.035224852 -0.000450088 
fit_wcls$beta_se_adjusted
# Intercept          Day 
# 0.0312540260 0.0006516683 
rbind(fit_wcls$beta_hat - t_quantile * fit_wcls$beta_se_adjusted,
      fit_wcls$beta_hat + t_quantile * fit_wcls$beta_se_adjusted)
# Intercept           Day
# [1,] -0.02784833 -0.0017652078
# [2,]  0.09829803  0.0008650318
2 * pt(abs(fit_wcls$beta_hat) / fit_wcls$beta_se_adjusted, 45 - 1 - 2, lower.tail = FALSE)
# Intercept       Day 
# 0.2661194 0.4935722 



### analysis 3: causal excursion effect moderation by gender


control_vars <- c("Day", "Y_previous_day", "gender")
moderator_vars <- "gender"

fit_wcls <- weighted_centered_least_square(
  dta = dta,
  id_varname = "study_id",
  decision_time_varname = "Day",
  treatment_varname = "message_sent",
  outcome_varname = "Y",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = 0.5,
  estimator_initial_value = NULL
)

t_quantile <- qt(0.975, 45 - 1 - 2) 

# estimator, SE, 95% CI, p-value
fit_wcls$beta_hat
# Intercept      gender 
# -0.00563512  0.02586462 
fit_wcls$beta_se_adjusted
# Intercept     gender 
# 0.01735056 0.03211170 
rbind(fit_wcls$beta_hat - t_quantile * fit_wcls$beta_se_adjusted,
      fit_wcls$beta_hat + t_quantile * fit_wcls$beta_se_adjusted)
# Intercept      gender
# [1,] -0.04064997 -0.03893942
# [2,]  0.02937973  0.09066866
2 * pt(abs(fit_wcls$beta_hat) / fit_wcls$beta_se_adjusted, 45 - 1 - 2, lower.tail = FALSE)
# Intercept    gender 
# 0.7469596 0.4250923 




### analysis 4: causal excursion effect moderation by Y_previous_day


control_vars <- c("Day", "Y_previous_day", "gender")
moderator_vars <- "Y_previous_day"

fit_wcls <- weighted_centered_least_square(
  dta = dta,
  id_varname = "study_id",
  decision_time_varname = "Day",
  treatment_varname = "message_sent",
  outcome_varname = "Y",
  control_varname = control_vars,
  moderator_varname = moderator_vars,
  rand_prob_varname = "prob_A",
  rand_prob_tilde_varname = NULL,
  rand_prob_tilde = 0.5,
  estimator_initial_value = NULL
)

t_quantile <- qt(0.975, 45 - 1 - 2) 

# estimator, SE, 95% CI, p-value
fit_wcls$beta_hat
# Intercept Y_previous_day 
# 0.017135845   -0.003479759 
fit_wcls$beta_se_adjusted
# Intercept Y_previous_day 
# 0.09521621     0.09406829 
rbind(fit_wcls$beta_hat - t_quantile * fit_wcls$beta_se_adjusted,
      fit_wcls$beta_hat + t_quantile * fit_wcls$beta_se_adjusted)
# Intercept Y_previous_day
# [1,] -0.1750182     -0.1933173
# [2,]  0.2092899      0.1863577
2 * pt(abs(fit_wcls$beta_hat) / fit_wcls$beta_se_adjusted, 45 - 1 - 2, lower.tail = FALSE)
# Intercept Y_previous_day 
# 0.8580434      0.9706668 