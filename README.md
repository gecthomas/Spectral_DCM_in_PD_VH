# Spectral DCM analysis in PD with VH (for SPM May 2024 practical session)
## This code contains a subset of the participants in the following paper: https://academic.oup.com/braincomms/article/5/1/fcac329/6895901
## It also contains simplified analyis and should NOT be expected to reproduce the results from the paper, it is for educational purposes only.
### DCM = dynamic causal model(ling)
### PD = Parkinson's disease 
### VH = visual hallucinations

In order to run the code in this directory, the dataset should be downloaded from the URL **[HERE](https://zenodo.org/records/11191180)** and unzipped into the master directory.

Repo contents (after unzipping)
  * **/code/**
    * **edit_SPM_paths.m** - Edits paths in the SPM.mat files to match local directory structure
    * **batch_extract_VOI_timeseries.m** - Extracts timeseries for 8 visual network nodes for all subjects
    * **run_first_level_DCM.m** - Specifies and estimates subject level spectral DCMs for all subjects
    * **run_PEB.m** - Parametric Empirical Bayes analysis to estimate group-level statistics from subject-level DCMs
    * **run_family_analysis.m** - Hypothesis based group-level analysis
    * **run_LOO.m** - Leave-one-out cross validation
    * **pointbiserial.m** - Code to test the point-biserial correlation *(credit to [Frederik Nagel](https://www.mathworks.com/matlabcentral/fileexchange/11222-point-biserial-correlation))*
  * **/func/**
    * **sub-001/**
      * **rsfmri_ica-aroma0000.nii**
      * ...
      * **rsfmri_ica-aroma0104.nii**
    * ...
    * **sub-053/**
      * **rsfmri_ica-aroma0000.nii**
      * ...
      * **rsfmri_ica-aroma0104.nii**
  * **/GLM/**
    * **sub-001/**
      * **beta_0001.nii**
      * **mask.nii**
      * **ResMS.nii**
      * **RPV.nii**
      * **SPM.mat**
    * ...
    * **sub-053/**
      * **beta_0001.nii**
      * **mask.nii**
      * **ResMS.nii**
      * **RPV.nii**
      * **SPM.mat**
  * **/analyses/**
    * **dmx_VH.mat** - design matrix containing subject sex, age, and VH status
    
