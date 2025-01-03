# Code to reproduce the figures from 'Analysis of multi-condition single-cell data with latent embedding multivariate regression'

> Ahlmann-Eltze C, Huber W (2025).
> “Analysis of multi-condition single-cell data with latent embedding multivariate regression.” Nature Genetics (2025).
> [doi:10.1038/s41588-024-01996-0](https://doi.org/10.1038/s41588-024-01996-0).

This repository contains the code to reproduce all figures from the manuscript. It is structured into 5 folders and several additional files which I will explain below.

-   `illustrations/`: contains the illustrator, png, and pdf files of the experimental design.
-   `plots/`: contains the finished figures from the main and supplementary text.
-   `notebooks/`: contains Rmarkdown files and rendered html notebooks which create the plots.
-   `benchmark/`: contains the compute intense code that was run on EMBL's cluster system for efficiency reasons.
-   `renv/activate.R`, `renv.lock`: the configuration files to recreate the computation environment in R (see [renv package](https://rstudio.github.io/renv/)).
-   `copy_benchmark_project.sh`: script that copies the files from the EMBL cluster system to my local computer.
-   `render_notebooks.sh`: script that calls rmarkdown's html rendering for each file in the `notebooks/` folder.

# Rendered Notebooks

-   [Preparation](https://htmlpreview.github.io/?https://github.com/const-ae/lemur-Paper/blob/master/benchmark/submission/prepare_glioblastoma.html) and [visualization](https://htmlpreview.github.io/?https://github.com/const-ae/lemur-Paper/blob/master/notebooks/glioblastoma_analysis.html) of panobinostat treatment in glioblastoma. 
-   [Preparation](https://htmlpreview.github.io/?https://github.com/const-ae/lemur-Paper/blob/master/benchmark/submission/prepare_alzheimer_plaques.html) and [visualization](https://htmlpreview.github.io/?https://github.com/const-ae/lemur-Paper/blob/master/notebooks/alzheimer_plaques_analysis.html) of the Alzheimer plaque densities
-   [Preparation](https://htmlpreview.github.io/?https://github.com/const-ae/lemur-Paper/blob/master/benchmark/submission/prepare_zebrafish.html) and [visualization](https://htmlpreview.github.io/?https://github.com/const-ae/lemur-Paper/blob/master/notebooks/zebrafish_analysis.html) of the zebrafish embryonic development
-   [Presentation](https://htmlpreview.github.io/?https://github.com/const-ae/lemur-Paper/blob/master/notebooks/benchmark_results.html) of the [benchmark](https://github.com/const-ae/lemur-Paper/blob/master/benchmark/submission/submit_data_preprocessing.R) results
-   [Visualization](https://htmlpreview.github.io/?https://github.com/const-ae/lemur-Paper/blob/master/notebooks/simulation_study.html) of a LEMUR fit on simulated data

# Software

The R package used in the analysis is available at https://bioconductor.org/packages/lemur/. The source code is at https://github.com/const-ae/lemur. 
A Python implementation of the LEMUR model is available at https://github.com/const-ae/pylemur.
