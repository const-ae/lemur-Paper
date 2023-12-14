
# LEMUR benchmark directory

This folder contains the scripts and folder structure for the benchmarking bits of the LEMUR paper. For efficiency reasons, all computationally intensive were run on our EMBL cluster system. The plotting was done locally though to iterate more quickly. The data in this folder copied over to my local machine by calling `../copy_benchmark_project.sh`.



# Tasks

1. Produce data for integration comparison (launched from `submit_data_preprocessing.R` with `MyWorkflowManager`)
2. Produce UMAP plots of all integrated datasets and save the Kang results to `"output/kang_visualization-{variation}.tsv"` (launched from `submit_data_preprocessing.R` with `MyWorkflowManager` and some custom script at the bottom of the script)
3. Produce the speed measurements and variance explained fraction
4. Prepare the biological data for local plotting:
  - Call `sbatch submission/run_prepare_glioblastoma.sh`
  - Call `sbatch submission/run_prepare_zebrafish.sh`
  - Call `sbatch submission/run_prepare_alzheimer_plaques.sh`
  




