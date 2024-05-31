## Description
This repository contains code to reproduce the analyses in Jones et al. (2024) "What makes a good tree metric?". We comapre a number of commonly used metrics on labelled phylogenetic trees. 

## Contents 
- Likelihood difference vs distance (Figures 4 and 5, and Supp figs): `likelihood-distance.R`
- MDS (Figure 3): `dist_mds_corr.R`
- Gisaid information for the alpha and delta data sets: `Cov2_GISAID_information.xlsx`

## Additional information 
We forked the RWTY R package to add additional metrics to the approximate and pseudo-ESS functions in RWTY. We have also added the ability to use any user-defined distance metric, where a
user-defined metric can be any function which takes two phylogenetic trees (from ape) as
input and outputs a number. 
Our changes are available as a fork of the original RWTY
repository at https://github.com/r8roy/RWTY
