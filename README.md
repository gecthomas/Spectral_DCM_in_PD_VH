# Spectral DCM analysis in PD with VH
## This code will reproduce the group-level results from: <placeholder>
### DCM = dynamic causal model(ling)
### PD = Parkinson's disease 
### VH = visual hallucinations


Repo contents
  * **/code/**
    * **run_PEB.m** - Parametric Empirical Bayes analysis to estimate group-level statistics from subject-level DCMs
    * **run_family_analysis.m** - Hypothesis based group-level analysis
    * **run_LOO.m** - Leave-one-out cross validation
    * **pointbiserial.m** - Code to test the point-biserial correlation *(credit to [Frederik Nagel](https://www.mathworks.com/matlabcentral/fileexchange/11222-point-biserial-correlation))*
  * **/analyses/**
    * **VN_LR_no_regr_GCM_full_VH.mat** - file containing estimated DCMs for all subjects
    * **dmx_VH.mat** - design matrix containing subject sex, age, and VH status
    
