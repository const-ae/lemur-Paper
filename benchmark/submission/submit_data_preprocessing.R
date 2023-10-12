setwd("/g/huber/users/ahlmanne/projects/lemur-Paper-benchmark")
library(tidyverse)
library(MyWorkflowManager)
source("submission/wrap_scripts.R")
# source("submission/job_submission_utils.R")
init("/scratch/ahlmanne/lemur_benchmark")

make_integration_jobs <- function(dataset, variation = "small"){
  job <- prepare_data(list(dataset = dataset,  variation = variation))  
  integration_jobs <- list(
    CPA = cpa_integration_prediction(list(data_id = job$result_id, dataset_config = dataset), dep_job = job),
    harmony = harmony_integration(list(data_id = job$result_id, dataset_config = dataset), dep_job = job),
    PCA = pca_integration_prediction(list(data_id = job$result_id, dataset_config = dataset), dep_job = job),
    lemur = lemur_integration_prediction(list(data_id = job$result_id, dataset_config = dataset), dep_job = job),
    invertible_harmony = lemur_integration_prediction(list(data_id = job$result_id, dataset_config = dataset, skip_multi_cond_pca = "true"), dep_job = job),
    multiCondPCA = lemur_integration_prediction(list(data_id = job$result_id, dataset_config = dataset, skip_alignment = "true"), dep_job = job)
  )
  eval_jobs <- map(integration_jobs, \(x){
    evaluate_integration(list(data_id = job$result_id, dataset_config = dataset, integration_id = x$result_id), dep_jobs = list(job, x))
  })
  eval_jobs
}

make_prediction_jobs <- function(dataset, variation = "small"){
  job <- prepare_data(list(dataset = dataset,  variation = variation))  
  prediction_jobs <- list(
    linear = linear_prediction(list(data_id = job$result_id, dataset_config = dataset), dep_job = job),
    CPA = cpa_integration_prediction(list(data_id = job$result_id, dataset_config = dataset), dep_job = job),
    PCA = pca_integration_prediction(list(data_id = job$result_id, dataset_config = dataset), dep_job = job),
    lemur = lemur_integration_prediction(list(data_id = job$result_id, dataset_config = dataset), dep_job = job),
    invertible_harmony = lemur_integration_prediction(list(data_id = job$result_id, dataset_config = dataset, skip_multi_cond_pca = "true"), dep_job = job),
    multiCondPCA = lemur_integration_prediction(list(data_id = job$result_id, dataset_config = dataset, skip_alignment = "true"), dep_job = job)
  )
  eval_jobs <- map(prediction_jobs, \(x){
    evaluate_prediction(list(data_id = job$result_id, dataset_config = dataset, prediction_id = x$result_id), dep_jobs = list(job, x))
  })
  eval_jobs
}


jobs1 <- make_prediction_jobs("kang", "random_holdout")
purrr::walk(jobs1, run_job, priority = "normal")
jobs2 <- make_integration_jobs("kang", "random_holdout")
purrr::walk(jobs2, run_job, priority = "normal")

map_chr(jobs2, job_status)
show_output_log(jobs1[["PCA"]]$dependencies[[1]])


all_jobs <- flatten(map(c("kang", "angelidis", "skinnider"), make_integration_jobs, "random_holdout"))

job_df <- tibble(name = names(all_jobs), job = all_jobs, 
                 data = rep(c("kang", "angelidis", "skinnider"), each = length(all_jobs)/3)) %>%
  mutate(status = map_chr(job, job_status))

tmp <- job_df %>%
  filter(status == "done") %>%
  mutate(res = map(job, \(j) MyWorkflowManager:::result_file_path(j) |> read_tsv(show_col_types = FALSE))) %>%
  mutate(stat = map(job, \(j) MyWorkflowManager:::stats_file_path(j) |> read_delim(show_col_types = FALSE, col_names = FALSE, delim = " ") |> deframe())) %>%
  unnest(res) %>%
  unnest_wider(stat, names_sep = "-")

write_tsv(tmp, "output/random_holdout_integration_results.tsv")
