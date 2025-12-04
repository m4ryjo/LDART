# R Code for Long-term Cognitive Effects of an Incremental Blood Pressure Intervention in a Mortal Cohort

This repository contains R code for:
- Analyses of the Betula data
- Simulation study described in:

**Josefsson, Karalija, and Daniels**  
https://doi.org/10.48550/arXiv.2311.00357

---

## Contents
- `Betula_analysis/` – Scripts for Betula data analyses
- `Simulation_study/` – Scripts for simulation study

---

## The R Package GcompBART
**GcompBART** is an R package for performing G-computation using Bayesian Additive Regression Trees.  
The package is available here:
- https://github.com/m4ryjo/GcompBART

GcompBART leverages the modeling capabilities of https://cran.r-project.org/package=SoftBart 
internally to fit flexible observed data models, with the optional use of sparsity-inducing 
priors for longitudinal data as described in:

> *Long-term cognitive effects of an incremental blood pressure intervention in a mortal cohort.*  
> Josefsson, Karalija, and Daniels  
> [https://doi.org/10.48550/arXiv.2311.00357](https://doi.org/10.48550/arXiv.2311.00357)

The package follows the same interface style as the SoftBart package for ease of use.

### Install GcompBART
To install GcompBART from GitHub:
```r
install.packages("remotes") # if not installed
remotes::install_github("m4ryjo/GcompBART")
```
## How to Run

### Option 1: From GitHub
1. Clone this repository:
```r
git clone https://github.com/m4ryjo/LDART.git
cd LDART
```

(You can replace LDART with a custom folder name if desired.)
2. Install required packages (including GcompBART).
3. Open R and set the working directory to the repository root.
4. Run scripts in:
- Betula_analysis/ for Betula data analyses
- Simulation_study/ for simulation study


### Option 2: From Journal Supplement

1. Download the ZIP file containing the code from the journal's supplementary materials.
2. Extract the files to a local directory.
3. Open R and set the working directory to the extracted folder.
4. Install required packages (including GcompBART):
install.packages("remotes")
remotes::install_github("m4ryjo/GcompBART")
5. Run scripts as described above.

## Notes
The Betula data are available upon request and approval by the project's stearing group. 
Missing data can be handled under both a MAR (Missing At Random) as well as a MNAR 
(Missing Not At Random) assumption. Deaths are truncated.
Blood pressure thresholds can easily be adjusted in the scripts (default = 140 mmHg).

## Citation
If you use this code, please cite:

> *Long-term cognitive effects of an incremental blood pressure intervention in a mortal cohort.*  
> Josefsson, Karalija, and Daniels  
> [https://doi.org/10.48550/arXiv.2311.00357](https://doi.org/10.48550/arXiv.2311.00357)


