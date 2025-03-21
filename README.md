# A R packge for Dose-finding Interval Designs

This repository contains a fully functional R package designed to implement various interval-based dose-finding designs, including 
- CCD (Continual Comparison Design)
- mTPI (Modified Toxicity Probability Interval)
- BOIN (Bayesian Optimal Interval Design)
- Keyboard Design (mTPI2)
- UMPBI (Uniformly Most Powerful Bayesian Interval Design)

This R package is based on the work of [Ruitao Lin](https://github.com/ruitaolin/IntervalDesign) with minor corrections, better documentation, and an additional approach to estimate posterior toxicity probability using Bayesian Logistic Regression Model (BLRM). 

## Overview
Phase I oncology trials aim to determine the maximum tolerated dose (MTD)—the dose level that induces dose-limiting toxicities (DLT) at a rate closest to the predefined target toxicity probability.

Interval designs have gained significant attention in phase I clinical trials due to their simplicity and strong finite-sample performance. These designs operate by comparing the observed toxicity rate (or the number of DLTs) against a predefined toxicity tolerance interval, guiding dose escalation and de-escalation decisions.

## Getting Started
Download the repo manually (as a .zip file) or clone it using Git.

```command
git clone https://github.com/ChuanhuiLiu/dosefinder
```
## Usage
Under the main directory, open the `dosefinder.Rproj` in an IDE such as RStudio as a package project.

Switch to development mode and load necessary packages by running the following commands in the R console
```rscript
library(usethis);library(devtools);library(ggplot2);library(rstanarm)
```

## Functions/Modules
`R` folder contain all the exported function. For any question, use `?function.name`:
* `get.rule`: Generate escalation and de-escalation rules and termination boundaries for several interval designs including the CCD, mTPI, BOIN, Keyboard and UMPBI designs.
* `estimate.prob`: Bayesian Logistic Regression Model (BLRM) for modeling the posterior toxicity probability
* `pava.transform`: Pool Adjacent Violator Algorithm (PAVA)
* `select.mtd`: Select Maximum Tolerated Dose (MTD) in a Dose-Finding Trial
* `simulate.summary`: Operating Characteristics for Bayesian Dose-Finding Designs via Simulations

## Example

For evaluating mTPI + PAVA method, run the following script
```rscript
p.true<-c(0.08,0.10,0.20,0.30,0.45,0.60)
summary.simulation(target=0.3,p.true=p.true,ncohort=12,cohortsize=3,design=2,type=1, cutoff.eli=0.95, ntrial = 1000)
```

Expected Output
```
$Selec_Percent
[1] "0.4"  "4.2"  "32.0" "46.1" "16.1" "1.1" 

$Avg_Patients
[1] "4.4"  "5.8"  "10.9" "10.0" "4.2"  "0.6" 

$Avg_Toxicities
[1] "0.3" "0.6" "2.2" "3.0" "1.9" "0.4"

$Summary
[1] "Total average toxicities: 8.4" "Total average patients: 36.0" 
```

To obtain a better result, as a main contribution of this project, use this for evaluating mTPI + BLRM method.
```rscript
p.true<-c(0.08,0.10,0.20,0.30,0.45,0.60)
summary.simulation(target=0.3,p.true=p.true,ncohort=12,cohortsize=3,design=2,type=2, cutoff.eli=0.95, ntrial = 20)
```
Expected Output
```rscript
$Selec_Percent
[1] "10.0" "0.0"  "25.0" "60.0" "5.0"  "0.0" 

$Avg_Patients
[1] "4.0"  "4.5"  "11.6" "13.8" "1.5"  "0.6" 

$Avg_Toxicities
[1] "0.3" "0.2" "1.9" "4.8" "0.8" "0.3"

$Summary
[1] "Total average toxicities: 8.3" "Total average patients: 36.0" 
```

## References
* Ruitao Lin's [Github Repo](https://github.com/ruitaolin/IntervalDesign)
* Guo W, Wang SJ, Yang S, et al. A Bayesian interval dose-finding design addressing Ockham’s razor: mTPI-2. Contemp Clin Trials 2017; 58: 23–33.
* Ivanova A, Flournoy N, Chung Y. Cumulative cohort design for dose finding. J Stat Plan Inference 2007; 137: 2316–2327.
* Ji Y, Liu P, Li Y, et al. A modified toxicity probability interval method for dose-finding trials. Clin Trials 2010; 7: 653–663.
* Lin R, Yin G. Uniformly most powerful Bayesian interval design for phase I dose-Finding trials.
* Liu S, Yuan Y. Bayesian optimal interval designs for phase I clinical trials. J R Stat Soc Ser C 2015; 64: 507–523.
* Yan F, Mandrekar SJ and Yuan Y. Keyboard: A Novel Bayesian toxicity probability interval design for phase I clinical trials. Clin Cancer Res 2017; 23: 3994–4003
