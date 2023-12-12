#!/bin/sh

#SBATCH --nodes=1 
#SBATCH --time=03:00:00
#SBATCH --mem=120GB
#SBATCH --job-name=render_prepare_glioblastoma
#SBATCH --output=tmp/slurm_logs/log-%x.%j.out
#SBATCH --output=tmp/slurm_logs/log-%x.%j.err
  
module load Pandoc/2.19
module load R-bundle-Bioconductor/3.16-foss-2022b-R-4.2.2
  
Rscript -e "rmarkdown::render('submission/prepare_glioblastoma.Rmd', output_format = 'html_document')"