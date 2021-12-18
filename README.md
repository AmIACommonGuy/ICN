# ICN
Interconnected Communities Network

The goal of this project is to apply the [Interconnected Communities Network](https://doi.org/10.1093/bioinformatics/btab047) (ICN) to cancer patient data from [The Cancer Genome Atlas](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) (TCGA) project in order to identify similar cancer types and potential subtypes. The writeup of our results can be found [here](). A minimum working example of the ICN framework and analysis pipeline can be found in the vignettes. Our methods are implemented into an R package `ICN`, which can be installed using the following command: 

## Installation Guide

```
library(devtools)
install_github("lwa19/ICN", build_vignettes = TRUE)
```

## Vignettes: Minimum Working Example

You may access all the vignettes and further analyses by running: 

```
browseVignettes('ICN')
```

## file structure

```
ICN
| R/                            - all of the written R functions and scripts
  | Auxiliary.R                 - all of the helper functions (inaccessible to users)
  | InterCon.R                  - R package function -- identification of significant edges
  | KLtest.R                    - R package function -- calculating and testing interconnectedness between communities using KL divergence
  | NICE_fast_cut.R             - R package function -- construct communities using K-means clustering
  | NICE_simulation.R           - R script -- simulate data and cluster using NICE algorithm
  | step1.R                     - R script -- execute NICE algorithm for input correlation matrix
  | step2.R                     - R script -- process NICE cluster output and perform KL divergence test between communities
  | step3.R                     - R script -- Detect interconnected edges within two communities
  | utils.R                     - R package function -- post processing function for NICE cluster output
| examples/                     - saved simulation data
  | corr.rds                    - simulated correlation matrix
  | nice.rds                    - NICE cluster output of simulated data
| vignettes/                    - R package vignettes
  | mwe.rmd                     - A minimum working example using simulation data
```


## Group Members
Ziming Huang (ziming@umich.edu)

Xin Luo (luosanj@umich.edu)

Yao Song (yaos@umich.edu)

Lijia Wang (lijiaw@umich.edu)

## Acknowledgements

We are grateful for the instruction and help from Professor Hui Jiang and the GSI Jingyi Zhai, and the materials taught in the course BIOSTAT 625. 
