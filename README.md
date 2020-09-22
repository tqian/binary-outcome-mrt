# Code for paper "Estimating time-varying causal excursion effect in mobile health with binary outcomes"

Code to reproduce results in the paper [Estimating time-varying causal excursion effect in mobile health with binary outcomes
](https://academic.oup.com/biomet/advance-article/doi/10.1093/biomet/asaa070/5901535) (Biometrika, 2020) by Tianchen Qian, Hyesun Yoo, Predrag Klasnja, Daniel Almirall, Susan A. Murphy

Tianchen Qian
2020.09.21


A detailed guide on how to reproduce the results in the paper is in Section 2. Section 1 explains the code structure (useful if you would like to use the code for your own purpose).

Please ignore [filename_mapping.txt](filename_mapping.txt). This file is used for my own reference, as the filenames are different on my local computer.


## 1. Code Structure

In file name alphabetical order:

* [barifit_analysis.R](barifit_analysis.R): Code to load, pre-process, and conduct analysis of the BariFit data. This code is used to generate results in Section 7 and Appendix H of the paper. One would need access to BariFit data including the following in order to run the code.
    * Bari-Fit data/barifit_csv_files/MRT_activity_suggestion_data.csv
    * Bari-Fit data/barifit_csv_files/barifit_demog_EMR.csv
    * Bari-Fit data/barifit_csv_files/barifit_weight_EMR.csv
    * Bari-Fit data/barifit_csv_files/food_track_analysis_data_correct.csv
    
* [dgm_binary_ar1_covariate.R](dgm_binary_ar1_covariate.R): a generative model for MRT data, where the time-varying covariate follows an AR(1) process. This function is *not* used in the simulations presented in the final paper.

* [dgm_binary_categorical_covariate.R](dgm_binary_categorical_covariate.R): a generative model for MRT data, where the time-varying covariate is categorical. This function is used in the simulations presented in the paper.

* [estimators_robust_adhocery.R](estimators_robust_adhocery.R): an early version of the EMEE estimator ("binary_outcome_moderated_effect()"). This function is *not* used in the final paper. Because the simulation code sources this file, I included this file here.

* [estimators.R](estimators.R): the main file that contains implementation of all estimators in the paper. A few notables ones are:
    * efficient_ee(): the original ECE estimator (Section 4 of the paper)
    * efficient_ee_twostep(): the ECE estimator with an alternative two-step implementation (Appendix G of the paper)
    * efficient_ee_modified_weight(): the ECE estimator with truncated weight (Section 6.1 of the paper, see Equation (12) for the form of weight truncation). This estimator is the ECE estimator used in the simulation in Section 6 of the paper.
    * weighted_centered_least_square(): the EMEE estimator (Section 5 of the paper)
    * The rest functions in [estimators.R](estimators.R) were used for exploratory purposes and are now obsolete.

* [simulation_consistency.R](simulation_consistency.R): code for reproducing the simulation results on consistency in Section 6.2 of the paper

* [simulation_efficiency_appendix.R](simulation_efficiency_appendix.R): code for reproducing the additional simulation results on efficiency in Appendix F of the paper

* [simulation_efficiency.R](simulation_efficiency.R): code for reproducing the simulation results on efficiency in Section 6.3 of the paper


## 2. Results in Paper and Corresponding Code

### Simulation on consistency (Section 6.2 in paper) 

Code is [simulation_consistency.R](simulation_consistency.R)

### Simulation on efficiency (Section 6.3 in paper)

Code is [simulation_efficiency.R](simulation_efficiency.R)

### Data analysis of BariFit (Section 7 in paper)

Code is [barifit_analysis.R](barifit_analysis.R). See code section "Section 7: analysis in main paper" therein.

This code is for illustration purpose only, as the BariFit data is not publicly available.

### Additional simulation on efficiency (Appendix F)

Code is [simulation_efficiency_appendix.R](simulation_efficiency_appendix.R)

### Alternative implementation of the ECE estimator (Appendix G)

The function that implements this estimator is "efficient_ee_twostep()" in [estimators.R](estimators.R).

The code to reproduce Table G.1 (comparing the two implementations of ECE) is not uploaded yet. (I'm looking for it in my local folders and will upload once I find it...)

### Exploratory plots of BariFit data (Appendix H)

Code is [barifit_analysis.R](barifit_analysis.R). See code section "Appendix H: exploratory plots" therein.
