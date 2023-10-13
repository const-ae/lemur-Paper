setwd("/g/huber/users/ahlmanne/projects/lemur-Paper-benchmark")
library(tidyverse)
library(MyWorkflowManager)
source("submission/wrap_scripts.R")
source("src/utils/config_helper.R")
init("/scratch/ahlmanne/lemur_benchmark")

datasets <- names(get_data_config())
names(datasets) <- datasets

make_integration_jobs <- function(dataset, variation = "small"){
  job <- prepare_data(list(dataset = dataset,  variation = variation)) 
  mem <- if(dataset %in% c("reyfman","mouse_gastrulation", "hrvatin"))  "80GB" else "40GB"
  
  default_params <- list(data_id = job$result_id, dataset_config = dataset)
  integration_jobs <- list(
    CPA = cpa_integration_prediction(default_params, dep_job = job, memory = mem),
    harmony = harmony_integration(default_params, dep_job = job, memory = mem),
    PCA = pca_integration_prediction(default_params, dep_job = job, memory = mem),
    lemur = lemur_integration_prediction(default_params, dep_job = job, memory = mem),
    invertible_harmony = lemur_integration_prediction(c(default_params, list(skip_multi_cond_pca = "true")), dep_job = job, memory = mem),
    multiCondPCA = lemur_integration_prediction(c(default_params, list(skip_alignment = "true")), dep_job = job, memory = mem)
  )
  eval_jobs <- map(integration_jobs, \(x){
    evaluate_integration(list(data_id = job$result_id, dataset_config = dataset, integration_id = x$result_id), dep_jobs = list(job, x))
  })
  eval_jobs
}

make_prediction_jobs <- function(dataset, variation = "small"){
  job <- prepare_data(list(dataset = dataset,  variation = variation))  
  mem <- if(dataset %in% c("reyfman","mouse_gastrulation", "hrvatin"))  "80GB" else "40GB"

  default_params <- list(data_id = job$result_id, dataset_config = dataset)
  prediction_jobs <- list(
    linear = linear_prediction(default_params, dep_job = job, memory = mem),
    CPA = cpa_integration_prediction(default_params, dep_job = job, memory = mem),
    PCA = pca_integration_prediction(default_params, dep_job = job, memory = mem),
    lemur = lemur_integration_prediction(default_params, dep_job = job, memory = mem),
    invertible_harmony = lemur_integration_prediction(c(default_params, list(skip_multi_cond_pca = "true")), dep_job = job, memory = mem),
    multiCondPCA = lemur_integration_prediction(c(default_params, list(skip_alignment = "true")), dep_job = job, memory = mem)
  )
  eval_jobs <- map(prediction_jobs, \(x){
    evaluate_prediction(list(data_id = job$result_id, dataset_config = dataset, prediction_id = x$result_id), dep_jobs = list(job, x))
  })
  eval_jobs
}



all_jobs <- list_flatten(map(datasets, make_prediction_jobs, "random_holdout"))
sub_jobs <- all_jobs[str_ends(names(all_jobs), "PCA")]
map_chr(all_jobs, job_status)

purrr::walk(all_jobs, run_job, priority = "normal")

job_df <- tibble(name = names(all_jobs), job = all_jobs) %>%
  separate(name, sep = "_", extra = "merge", into = c("data", "method")) %>%
  mutate(status = map_chr(job, job_status))

tmp <- job_df %>%
  filter(status == "done") %>%
  mutate(res = map(job, \(j) MyWorkflowManager:::result_file_path(j) |> read_tsv(show_col_types = FALSE))) %>%
  mutate(stat = map(job, \(j) MyWorkflowManager:::stats_file_path(j) |> read_delim(show_col_types = FALSE, col_names = FALSE, delim = " ") |> deframe())) %>%
  unnest(res) %>%
  unnest_wider(stat, names_sep = "-")

write_tsv(tmp, "output/random_holdout_integration_results.tsv")
