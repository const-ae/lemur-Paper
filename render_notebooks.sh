#!/bin/sh

Rscript -e "rmarkdown::render('notebooks/benchmark_results.Rmd', output_format = 'html_document')"
Rscript -e "rmarkdown::render('notebooks/alzheimer_plaques_analysis.Rmd', output_format = 'html_document')"
Rscript -e "rmarkdown::render('notebooks/zebrafish_analysis.Rmd', output_format = 'html_document')"
Rscript -e "rmarkdown::render('notebooks/glioblastoma_analysis.Rmd', output_format = 'html_document')"
Rscript -e "rmarkdown::render('notebooks/simulation_study.Rmd', output_format = 'html_document')"
