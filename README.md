# External Validation of EPIC-US using NOVELTY

## Overview
This repository contains the codebase for the external validation of the EPIC-US to the the NOVELTY observational study. 

## Repository Structure & Key Files

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
* `epicUS` package
* Standard data manipulation and visualization libraries (e.g. `tidyverse`, `ggplot2`, `scales`, `ggthemes`)

## Usage
Ensure the `epicUS` package is properly installed prior to running the 
calibration scripts and vignettes.

## Installation Process
### Windows 7 or Later
1. Download and Install the latest version of R from [https://cran.r-project.org/bin/windows/base/](https://cran.r-project.org/bin/windows/base/)
2. Download and Install R Studio from [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/)
3. Download and Install the latest version of Rtools from [https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/) 
4. Using either an R session in Terminal or in R Studio, install the package `devtools`:

```r
  install.packages('remotes')
```

5. Install epicR from GitHub:

```r
remotes::install_github('resplab/epicUS')
```

### Mac OS Sierra and Later
1. Download and Install the latest version of R from [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/)
2. Download and Install R Studio from [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/)
3. Open the Terminal and remove previous installations of `clang`:

```bash
# Delete the clang4 binary
sudo rm -rf /usr/local/clang4
# Delete the clang6 binary
sudo rm -rf /usr/local/clang6

# Delete the prior version of gfortran installed
sudo rm -rf /usr/local/gfortran
sudo rm -rf /usr/local/bin/gfortran

# Remove the install receipts that indicate a package is present

# Remove the gfortran install receipts (run after the above commands)
sudo rm /private/var/db/receipts/com.gnu.gfortran.bom
sudo rm /private/var/db/receipts/com.gnu.gfortran.plist

# Remove the clang4 installer receipt
sudo rm /private/var/db/receipts/com.rbinaries.clang4.bom
sudo rm /private/var/db/receipts/com.rbinaries.clang4.plist

# Remove the Makevars file
rm ~/.R/Makevars
```
4. Install the latest version of `clang` by installing Xcode command tools: 
`xcode-select --install

5. Install the appropriate version of `gfortran` based on your Mac OS version using the dmg file found at [https://github.com/fxcoudert/gfortran-for-macOS/releases](https://github.com/fxcoudert/gfortran-for-macOS/releases) 

6. Using either an R session in Terminal or in R Studio, install the packages `remotes` and `usethis`:

```r
install.packages (c('remotes', 'usethis'))
```
7. Open the `.Renviron` file usibng the command `usethis::edit_r_environ()`
8. Add the line `PATH="/usr/local/clang8/bin:${PATH}"` to the file. If you installed any clang version above 8, modify the file accordingly. Save the `.Renviron` file and restart R.  

9. You should now be able to Install epicR from GitHub:
```r
remotes::install_github('resplab/epicR')
```

### Ubuntu 22.04 and Later
1. Install R by executing the following commands in Terminal:

```bash
  sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
```
```bash
  sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
```
```bash
  sudo apt-get update
```
```bash
  sudo apt-get install r-base
```
If the installation is successful, you should be able to start R:
```bash
  sudo -i R
```

2. Download and Install R Studio from [https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/)
3. Install `libcurl` from Terminal: 

```bash
  sudo apt-get install libcurl4-openssl-dev libssl-dev r-base-dev
```

4. Using either an R session in Terminal or in R Studio, install the package `devtools`:

```r
install.packages ('remotes')
```
  
5. Install epicR from GitHub:

```r
remotes::install_github('resplab/epicR')
```

# Quick Guide

To run EPIC with default inputs and settings, use the code snippet below. 
```
library(epicR)
init_session()
run()
Cget_output()
terminate_session()
```
Default inputs can be retrieved with `get_input()`, changed as needed, and resubmitted as a parameter to the run function:
```
init_session()
input <- get_input()
input$values$global_parameters$time_horizon <- 5
run(input=input$values)
results <- Cget_output()
resultsExra <- Cget_output_ex()
terminate_session()

```

For some studies, having access to the entire event history of the simulated population might be beneficial. Capturing event history is possible by setting  `record_mode` as a `setting`. 

```
settings <- get_default_settings()
settings$record_mode <- 2
settings$n_base_agents <- 1e4
init_session(settings = settings)
run()
results <- Cget_output()
events <- as.data.frame(Cget_all_events_matrix())
head(events)
terminate_session()

```
Note that you might need a large amount of memory available, if you want to collect event history for a large number of patients. 

In the events data frame, each type of event has a code corresponding to the table below:

|Event|No.|
|-----|---|
|start |0 |
|annual|1 |
|birthday| 2 |
|smoking change | 3|
|COPD incidence | 4|
|Exacerbation | 5 |
|Exacerbation end| 6|
|Death by Exacerbation | 7|
|Doctor visit | 8|
|Medication change | 9|
|Background death | 13|
|End | 14|
