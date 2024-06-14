setwd("/g/huber/users/ahlmanne/projects/lemur-Paper-benchmark")
library(tidyverse)
library(MyWorkflowManager)
source("submission/wrap_scripts.R")
source("src/utils/config_helper.R")
init("/scratch/ahlmanne/lemur_benchmark")

datasets <- names(get_data_config())
names(datasets) <- datasets

make_integration_jobs <- function(dataset, variation = "small"){
  mem <- if(dataset %in% c("reyfman","mouse_gastrulation", "hrvatin", "goldfarbmuren"))  "160GB" else "40GB"
  data_job <- prepare_data(list(dataset = dataset,  variation = variation), memory = mem) 
  
  default_params <- list(data_id = data_job$result_id, dataset_config = dataset)
  ref_jobs <- list(
    CPA_ref = cpa_on_counts_integration_prediction(c(default_params, list("no_integration"="")), dep_job = data_job, memory = mem),
    scvi_ref = scvi_integration_prediction(c(default_params, list("no_integration"="")), dep_job = data_job, memory = mem),
    PCA_ref = pca_integration_prediction(default_params, dep_job = data_job, memory = mem)
  )
  
  integration_jobs <- list(
    CPA_kangparams = cpa_on_counts_integration_prediction(default_params, dep_job = data_job, memory = mem),
    scvi = scvi_integration_prediction(default_params, dep_job = data_job, memory = mem),
    harmony = harmony_integration(default_params, dep_job = data_job, memory = mem),
    PCA = pca_integration_prediction(default_params, dep_job = data_job, memory = mem),
    lemur = lemur_integration_prediction(default_params, dep_job = data_job, memory = mem),
    invertible_harmony = lemur_integration_prediction(c(default_params, list(skip_multi_cond_pca = "true")), dep_job = data_job, memory = mem),
    multiCondPCA = lemur_integration_prediction(c(default_params, list(skip_alignment = "true")), dep_job = data_job, memory = mem)
  )
  eval_jobs <- map2(integration_jobs, names(integration_jobs), \(job, n){
    ref_job <- if(n == "CPA_kangparams") ref_jobs[["CPA_ref"]]
    else if(n == "scvi") ref_jobs[["scvi_ref"]]
    else ref_jobs[["PCA_ref"]]
    evaluate_integration(list(data_id = data_job$result_id, dataset_config = dataset, 
                              integration_id = job$result_id, reference_id = ref_job$result_id),
                         dep_jobs = list(data_job, job, ref_job))
  })
  eval_jobs
}

make_visualization_jobs <- function(dataset, variation = "small"){
  mem <- if(dataset %in% c("reyfman","mouse_gastrulation", "hrvatin", "goldfarbmuren"))  "160GB" else "40GB"
  job <- prepare_data(list(dataset = dataset,  variation = variation), memory = mem) 
  
  default_params <- list(data_id = job$result_id, dataset_config = dataset)
  integration_jobs <- list(
    CPA_kangparams = cpa_on_counts_integration_prediction(default_params, dep_job = job, memory = mem),
    scvi = scvi_integration_prediction(default_params, dep_job = job, memory = mem),
    harmony = harmony_integration(default_params, dep_job = job, memory = mem),
    PCA = pca_integration_prediction(default_params, dep_job = job, memory = mem),
    lemur = lemur_integration_prediction(default_params, dep_job = job, memory = mem),
    invertible_harmony = lemur_integration_prediction(c(default_params, list(skip_multi_cond_pca = "true")), dep_job = job, memory = mem),
    multiCondPCA = lemur_integration_prediction(c(default_params, list(skip_alignment = "true")), dep_job = job, memory = mem)
  )
  vis_jobs <- map(integration_jobs, \(x){
    visualize_integration(list(data_id = job$result_id, dataset_config = dataset, integration_id = x$result_id), dep_jobs = list(job, x))
  })
  vis_jobs
}

make_prediction_jobs <- function(dataset, variation = "small"){
  mem <- if(dataset %in% c("reyfman","mouse_gastrulation", "hrvatin", "goldfarbmuren"))  "160GB" else "40GB"
  job <- prepare_data(list(dataset = dataset,  variation = variation), memory = mem) 

  default_params <- list(data_id = job$result_id, dataset_config = dataset)
  prediction_jobs <- list(
    linear = linear_prediction(default_params, dep_job = job, memory = mem),
    no_change = no_change_prediction(default_params, dep_job = job, memory = mem),
    CPA_large = cpa_integration_prediction(c(default_params, list(n_latent = 64, max_epochs = 2000)), dep_job = job, memory = mem),
    CPA_kangparams = cpa_kang_params_integration_prediction(default_params, dep_job = job, memory = mem),
    scvi = scvi_integration_prediction(default_params, dep_job = job, memory = mem),
    PCA = pca_integration_prediction(default_params, dep_job = job, memory = mem),
    lemur = lemur_integration_prediction(default_params, dep_job = job, memory = mem),
    invertible_harmony = lemur_integration_prediction(c(default_params, list(skip_multi_cond_pca = "true")), dep_job = job, memory = mem),
    multiCondPCA = lemur_integration_prediction(c(default_params, list(skip_alignment = "true")), dep_job = job, memory = mem)
  )
  eval_jobs <- map(prediction_jobs, \(x){
    evaluate_prediction(list(data_id = job$result_id, dataset_config = dataset, prediction_id = x$result_id), dep_jobs = list(job, x),
                        memory = mem)
  })
  eval_jobs
}

make_variance_explained_jobs <- function(dataset, variation = "small"){
  mem <- if(dataset %in% c("reyfman","mouse_gastrulation", "hrvatin", "goldfarbmuren"))  "120GB" else "40GB"
  job <- prepare_data(list(dataset = dataset,  variation = variation), memory = mem) 
  pca_dims <- c(1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 50, 75, 100, 150, 200)
  calc_variance_explained(params = list(data_id = job$result_id, dataset_config = dataset, pca_dims = pca_dims),
                          dep_jobs = list(job), memory = mem)
}

variation <- "random_holdout_hvg"
data_jobs <- map(datasets, \(dat) prepare_data(list(dataset = dat,  variation = variation), memory = "40GB"))
map_chr(data_jobs, job_status)
walk(data_jobs, run_job, priority = "normal")

int_jobs <- list_flatten(map(datasets, make_integration_jobs, variation), name_spec = "{outer}-{inner}")
stat <- map_chr(int_jobs, job_status); table(stat)
write_rds(int_jobs, glue::glue("tmp/integration_jobs-{variation}.RDS"))
walk(int_jobs, run_job, priority = "normal")

pred_jobs <- list_flatten(map(datasets, make_prediction_jobs, variation), name_spec = "{outer}-{inner}")
stat <- map_chr(pred_jobs, job_status); table(stat)
write_rds(pred_jobs, glue::glue("tmp/prediction_jobs-{variation}.RDS"))
walk(pred_jobs, run_job, priority = "normal")

vis_jobs <- list_flatten(map(datasets, make_visualization_jobs, variation), name_spec = "{outer}-{inner}")
stat <- map_chr(vis_jobs, job_status); table(stat)
write_rds(vis_jobs, glue::glue("tmp/visualization_jobs-{variation}.RDS"))
walk(vis_jobs, run_job, priority = "normal")

var_jobs <- c(list_flatten(map(datasets, make_variance_explained_jobs, "random_holdout_hvg"), name_spec = "{outer}-{inner}"),
              list_flatten(map(datasets, make_variance_explained_jobs, "random_holdout"), name_spec = "{outer}-{inner}"))
stat <- map_chr(var_jobs, job_status); table(stat)
write_rds(var_jobs, glue::glue("tmp/variance_explained_jobs.RDS"))
walk(var_jobs, run_job, priority = "normal")
var_jobs %>%
  keep(\(j) job_status(j) == "done") %>%
  map_df(\(j) read_tsv(result_file_path(j), show_col_types = FALSE)) %>%
  write_tsv(glue::glue("output/variance_explained.tsv"))

get_result_table <- function(jobs){
  job_df <- tibble(name = names(jobs), job = jobs) %>%
    separate(name, sep = "-", into = c("data", "method")) %>%
    mutate(status = map_chr(job, job_status))
  
  job_df %>%
    filter(status == "done") %>%
    mutate(res = map(job, \(j) MyWorkflowManager:::result_file_path(j) |> read_tsv(show_col_types = FALSE))) %>%
    mutate(stat = map(job, \(j) MyWorkflowManager:::stats_file_path(j) |> read_delim(show_col_types = FALSE, col_names = FALSE, delim = " ") |> deframe())) %>%
    unnest(res) %>%
    unnest_wider(stat, names_sep = "-")
}

int_tab <- get_result_table(int_jobs)
write_tsv(int_tab, glue::glue("output/integration_results-{variation}.tsv"))

pred_tab <- get_result_table(pred_jobs)
write_tsv(pred_tab, glue::glue("output/prediction_results-{variation}.tsv"))


# -------------
# Make plots
vis_df <- tibble(name = names(vis_jobs), job = vis_jobs) %>%
  separate(name, sep = "-", into = c("data", "method")) %>%
  mutate(status = map_chr(job, job_status)) %>%
  filter(status == "done") %>%
  mutate(vis = map(job, \(j) read_tsv(MyWorkflowManager:::result_file_path(j), show_col_types = FALSE)))

plots <- vis_df %>%
  # filter(data %in% c("angelidis", "aztekin", "kang")) %>%
  # filter(as.integer(as.factor(data)) <  12 ) %>%
  group_by(data) %>%
  group_map(\(dat, key){
    dat %>%
      unnest(vis) %>%
      # sample_n(size = min(nrow(dat$vis[[1]]), 1e4)) %>%
      ggplot(aes(x = umap1, y = umap2)) +
      ggrastr::rasterise(geom_point(aes(color = covariate), size = 0.3, stroke = 0), dpi = 300) +
      # geom_point(aes(color = covariate), size = 0.1, stroke = 0) +
      facet_wrap(vars(method), nrow = 1) +
      coord_fixed() +
      labs(title = key[[1]])
  })
pl <- cowplot::plot_grid(plotlist = plots, nrow = length(plots))
cowplot::save_plot(glue::glue("tmp/embedding_vis-{variation}.pdf"), pl, ncol = 6, nrow = length(plots), base_asp = 1)

tibble(name = names(vis_jobs), job = vis_jobs) %>%
  separate(name, sep = "-", into = c("data", "method")) %>%
  filter(data == "kang") %>%
  mutate(status = map_chr(job, job_status)) %>%
  filter(status == "done") %>%
  mutate(vis = map(job, \(j) read_tsv(MyWorkflowManager:::result_file_path(j), show_col_types = FALSE))) %>%
  unnest(vis) %>%
  write_tsv(glue::glue("output/kang_visualization-{variation}.tsv"))

# -------------
# Kang predictions
library(MatrixGenerics)
library(SingleCellExperiment)
kang_pred <- make_prediction_jobs("kang", variation = variation)
map_chr(kang_pred, job_status)
holdout_h5ad_file <- file.path(result_file_path(kang_pred$lemur$dependencies[[1]]), "holdout.h5ad")
holdout_sce <- zellkonverter::readH5AD(holdout_h5ad_file)
is_stim <- holdout_sce$label == "stim"

celltype_split <- lapply(unique(holdout_sce$cell_type), \(ct) which(holdout_sce$cell_type == ct))
names(celltype_split) <- unique(holdout_sce$cell_type)

cell_type_pred <- bind_rows(lapply(names(kang_pred), \(method){
  holdout_file1 <- file.path(result_file_path(kang_pred[[method]]$dependencies[[2]]), "holdout-prediction_ctrl.tsv")
  pred1 <- as.matrix(data.table::fread(holdout_file1, sep = "\t", header = FALSE))
  holdout_file2 <- file.path(result_file_path(kang_pred[[method]]$dependencies[[2]]), "holdout-prediction_stim.tsv")
  pred2 <- as.matrix(data.table::fread(holdout_file2, sep = "\t", header = FALSE))
  pred_stim <- lemur:::aggregate_matrix(pred2, group_split = celltype_split, aggr_fnc = rowMeans2, col_sel = !is_stim)
  pred_ctrl <- lemur:::aggregate_matrix(pred1, group_split = celltype_split, aggr_fnc = rowMeans2, col_sel = is_stim)
  colnames(pred_stim) <- paste0("stim-", colnames(pred_stim))
  colnames(pred_ctrl) <- paste0("ctrl-", colnames(pred_ctrl))
  
  tibble(method = method, gene = rownames(holdout_sce)) %>%
    bind_cols(as_tibble(pred_stim), as_tibble(pred_ctrl))
}))

obs_stim <- lemur:::aggregate_matrix(logcounts(holdout_sce), group_split = celltype_split, aggr_fnc = rowMeans2, col_sel = is_stim)
obs_ctrl <- lemur:::aggregate_matrix(logcounts(holdout_sce), group_split = celltype_split, aggr_fnc = rowMeans2, col_sel = !is_stim)
colnames(obs_stim) <- paste0("stim-", colnames(obs_stim))
colnames(obs_ctrl) <- paste0("ctrl-", colnames(obs_ctrl))
cell_type_pred <- bind_rows(cell_type_pred,
                            tibble(method = "obs",gene = rownames(holdout_sce)) %>%  
                              bind_cols(as_tibble(obs_stim), as_tibble(obs_ctrl)))

write_tsv(cell_type_pred, glue::glue("output/kang-detailed_prediction_results_{variation}.tsv"))


