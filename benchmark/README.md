
# LEMUR benchmark directory

This folder contains the scripts and folder structure for the benchmarking bits of the LEMUR paper. For efficiency reasons, all computationally intensive were run on our EMBL cluster system. The plotting was done locally though to iterate more quickly. The data in this folder copied over to my local machine by calling `../copy_benchmark_project.sh`.

# Execution

The project is designed to execute all scripts using the [`MyWorkflowManager`](https://github.com/const-ae/MyWorkflowManager) and the slurm cluster manager. Install the `MyWorkflowManager` R package from Github and then follow the steps outlined in Tasks.
Alternatively, you can call each script directly from the command line providing the appropriate input arguments. For example to prepare a dataset for the differential expression benchmark, you could call:

```bash
ahlmanne@cluster$ Rscript src/prepare_data_for_de.R  \
  --dataset angelidis \
  --condition 3m \
  --randomize cells \
  --cut_at 1 10 25 \
  --n_de_genes 200 \
  --working_dir /tmp/ahlmanne/benchmark \
  --result_id 1
```

To subsequently run LEMUR's differential expression analysis, call:
```bash
ahlmanne@cluster$ Rscript src/lemur_differential_expression.R  \
  --data_id 1 \
  --dataset_config angelidis \
  --split_test_training FALSE \
  --working_dir /tmp/ahlmanne/benchmark \
  --result_id 1-split_test_training
```

# Tasks

1. Produce data for integration comparison (launched from `submit_data_preprocessing.R` with `MyWorkflowManager`)
2. Produce UMAP plots of all integrated datasets and save the Kang results to `"output/kang_visualization-{variation}.tsv"` (launched from `submit_data_preprocessing.R` with `MyWorkflowManager` and some custom script at the bottom of the script)
3. Produce the speed measurements and variance explained fraction
4. Prepare the biological data for local plotting:
  - Call `sbatch submission/run_prepare_glioblastoma.sh`
  - Call `sbatch submission/run_prepare_zebrafish.sh`
  - Call `sbatch submission/run_prepare_alzheimer_plaques.sh`
  




