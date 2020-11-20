# Optimal timing for cancer screening and adaptive surveillance using mathematical modeling

These are the scripts written and used to implement analyses and calculate solutions for the optimal timing framework, with applications using the multistage clonal espansion for esophageal adenocarcinoma (MSCE-EAC) model, that were used as part of the manuscript:

"Optimal timing for cancer screening and adaptive surveillance using mathematical modeling"

Authors:  Kit Curtius, Anup Dewanji, William Hazelton, Joel H Rubenstein, E. Georg Luebeck



For correspondence:
- Kit Curtius

Email: kcurtius@health.ucsd.edu


_Data used in the analyses are publicly available and previously published. Clinical Outcomes Research Initiative values from the paper were published in Rubenstein et al. 2010 Ref. 37_


### BE_density_Optimal_Timing.R
- Functions to compute density and cumulative risk functions for Barrett's esophagus (precursor to EAC) onset for given age

### am_Bootstrap_sensitivity_2020.R
- Parameters for stochastic rates of normal to EAC development in men drawn from MCMC posterior estimates (95% CI and SD computed)

### af_Bootstrap_sensitivity_2020.R
- Parameters for stochastic rates of normal to EAC development in women drawn from MCMC posterior estimates (95% CI and SD computed)

### Optimal_timing_screening.R
-R script to reproduce model predictions provided in manuscript Results including:
## Optimal screening ages for BE in men and women, general population and GERD-specific
## Optimal re-screening ages for BE in men and women, general population and GERD-specific
## Cancer risk associated with surveillance of BE in men and women
-Also provides code to reproduce model-predictions for Figures 2-5, Table 1, and Figure S1B
 
-Some scripts are purposed from other Github repositories, as indcated in the file names.
