# Model-robust and efficient covariate adjustment for cluster-randomized experiments

## Data

### Abstract 
Three completed cluster-randomized experiments were analyzed in the paper for demonstration. The three data sets are from the Work, Family, and Health Study (WFHS), Pain Program for Active Coping and Training study (PPACT), and the Improving rational use of artemisinin combination therapies through diagnosis-dependent subsidies (ACTS) study.


### Availability 
The WFHS data are publicly available at https://www.icpsr.umich.edu/web/DSDR/studies/36158. The ACTS data are publicly available at https://datadryad.org/stash/dataset/doi:10.5061/dryad.59p4111. However, the PPACT data are not publicly available. Interested readers can contact the study team of paper "A Primary Care–Based Cognitive Behavioral Therapy Intervention for Long-Term Opioid Users With Chronic Pain" for further information.


## Code

### Abstract
The R code here can be used to reproduce all of simulations and data applications (if the data are ready).

### Description 
The folder `simulations` contains R code for our simulations.

 - `simulation-1-1.R` contains the code to reproduce Table 1 of the main paper. 
 - `simulation-1-2.R` contains the code to reproduce Table S1 of the main paper. 
 - `simulation-2-1.R` contains the code to reproduce Table 2 of the main paper. 
 - `simulation-2-2.R` contains the code to reproduce Table S2 of the main paper. 
 - `summary-simulation.R` contains the code to summarize all simulation results.
 
 The folder `data-analysis` contains R code for our data applications
 
 - `data-analysis-WFHS.R` contains the code for re-analyzing the WFHS data set.
 - `data-analysis-PPACT.R` contains the code for re-analyzing the PPACT data set.
 - `data-analysis-ACTS.R` contains the code for re-analyzing the ACTS data set.
 - `WHFS-resampled-analysis.R` contains the code for simulations based on the WFHS data set.

### Reproducibility 
Tables and Figures in the main text can be reproduced using the R scripts.



