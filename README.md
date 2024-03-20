# Code to reproduce the figures from 'Analysis of multi-condition single-cell data with latent embedding multivariate regression'

> Citation: Ahlmann-Eltze and Huber (2024), Analysis of multi-condition single-cell data with latent embedding multivariate regression.” bioRxiv. [doi:10.1101/2023.03.06.531268](http://dx.doi.org/10.1101/2023.03.06.531268)

This repository contains the code to reproduce all figures from the manuscript. It is structured into 5 folders and several additional files which I will explain below.

-   `illustrations/`: contains the illustrator, png, and pdf files of the experimental design.
-   `plots/`: contains the finished figures from the main and supplementary text.
-   `notebooks/`: contains Rmarkdown files and rendered html notebooks which create the plots.
-   `benchmark/`: contains the compute intense code that was run on EMBL's cluster system for efficiency reasons.
-   `renv/activate.R`, `renv.lock`: the configuration files to recreate the computation environment in R (see [renv package](https://rstudio.github.io/renv/)).
-   `copy_benchmark_project.sh`: script that copies the files from the EMBL cluster system to my local computer.
-   `render_notebooks.sh`: script that calls rmarkdown's html rendering for each file in the `notebooks/` folder.

# Rendered Notebooks

-   [Analysis of panobinostat treatment in glioblastoma](https://htmlpreview.github.io/?https://github.com/const-ae/lemur-Paper/blob/master/notebooks/glioblastoma_analysis.html)
-   [Analysis of the Alzheimer plaque densities](https://htmlpreview.github.io/?https://github.com/const-ae/lemur-Paper/blob/master/notebooks/alzheimer_plaques_analysis.html)
-   [Analysis of the zebrafish embryonic development](https://htmlpreview.github.io/?https://github.com/const-ae/lemur-Paper/blob/master/notebooks/zebrafish_analysis.html)
-   [Presentation of the benchmark results](https://htmlpreview.github.io/?https://github.com/const-ae/lemur-Paper/blob/master/notebooks/benchmark_results.html)
-   [LEMUR fit on simulated data](https://htmlpreview.github.io/?https://github.com/const-ae/lemur-Paper/blob/master/notebooks/simulation_study.html)

# Software

The R package used in the analysis is available at https://bioconductor.org/packages/lemur/. The source code is at https://github.com/const-ae/lemur. 
A Python implementation of the LEMUR model is available at https://github.com/const-ae/pylemur.
