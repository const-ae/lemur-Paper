setwd("/g/huber/users/ahlmanne/projects/lemur-Paper-benchmark")
library(tidyverse)
library(MyWorkflowManager)
source("submission/wrap_scripts.R")
# source("submission/job_submission_utils.R")
init("/scratch/ahlmanne/lemur_benchmark")


job <- prepare_data(list(dataset = "kang",  variation = "small"))
run_job(job, "normal")
job_status(job)

job2 <- linear_prediction(list(data_id = job$result_id, dataset_config = "kang"), dep_job = job)
run_job(job2, "normal")

py_job <- cpa_integration_prediction(list(data_id = job$result_id, dataset_config = "kang"), dep_job = job)
run_job(py_job, priority = "high")

harm_job <- harmony_integration(list(data_id = job$result_id, dataset_config = "kang"), dep_job = job)
run_job(harm_job, priority = "high")

pca_job <- pca_integration_prediction(list(data_id = job$result_id, dataset_config = "kang"), dep_job = job)
run_job(pca_job, priority = "high")
job_status(pca_job)
show_output_log(pca_job)
