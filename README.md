# External Validation of EPIC-US to the NOVELTY Trial

## Overview
This repository contains the codebase for the external validation of the EPIC-US to the the NOVELTY observational study. 

## Repository Structure and Key Files

* **`Copula Model.R`**
  This details the joint distribution of patient characteristics in NOVELTY 
  and creates a matched cohort of EPIC-US simulated patients. 
   
* **`NOVELTY vs. EPIC Matching and Results_LLN.md`**
  This details the codebase and the cohort matching and
  exacerbation model results between EPIC-US and NOVELTY. 
  
* **`RDS files`**
  These are object files:
  novelty_COPDLLN.rds describes the joint distribution of patient characteristics in NOVELTY
  epic_selected.rds and combined_histories.rds are the simulated EPIC-US patients used in analysis

## Requirements
* [epicUS](https://github.com/resplab/epicUS) package
* Standard data manipulation and visualization libraries (e.g. `tidyverse`, `ggplot2`, `scales`, `ggthemes`)

## Usage
Ensure the [epicUS](https://github.com/resplab/epicUS) package is properly installed prior to running the 
calibration scripts and vignettes.
