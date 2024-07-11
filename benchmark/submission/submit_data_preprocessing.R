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

make_perturbation_jobs <- function(datasets = c('norman', 'adamson', 'replogle_k562_essential', 'replogle_rpe1_essential'),
                                   seeds = 1:5){
  jobs <- lapply(datasets, \(dataset){
    inner_jobs <- lapply(seeds, \(se){
      config_job <- prepare_perturbation_data(list(dataset = dataset,  seed = se))
      default_params <- list(dataset_name = dataset, test_train_config_id = config_job$result_id, seed = se)
      pert_jobs <- list(
        scgpt = scgpt_combinatorial_prediction(default_params, dep_jobs = list(config_job)),
        scgpt_oneshot = scgpt_combinatorial_prediction(c(default_params, list(epochs = 0)), dep_jobs = list(config_job)),
        gears = gears_combinatorial_prediction(default_params, dep_jobs = list(config_job)),
        ground_truth = ground_truth_combinatorial_prediction(default_params, dep_jobs = list(config_job))
      )
      if(dataset == "norman"){
        pert_jobs <- append(pert_jobs, list(
          additive_model = additive_model_combinatorial_prediction(default_params, dep_jobs = list(config_job)),
          pylemur = pylemur_combinatorial_prediction(default_params, dep_jobs = list(config_job))
        ))
      }else{
        pert_jobs <- append(pert_jobs, list(
          linear = linear_perturbation_prediction(c(default_params, list(ridge_penalty = 0.1, pca_dim = 5)), dep_jobs = list(config_job)),
          linear_highRidge = linear_perturbation_prediction(c(default_params, list(ridge_penalty = 3, pca_dim = 5)), dep_jobs = list(config_job)),
          pretrained_k562 = transfer_perturbation_prediction(c(default_params, list(ridge_penalty = 0.1, pca_dim = 5, reference_data = "replogle_k562_essential")),  dep_jobs = list(config_job)),
          pretrained_k562_highDim = transfer_perturbation_prediction(c(default_params, list(ridge_penalty = 0.1, pca_dim = 40, reference_data = "replogle_k562_essential")),  dep_jobs = list(config_job)),
          pretrained_rpe1 = transfer_perturbation_prediction(c(default_params, list(ridge_penalty = 0.1, pca_dim = 5, reference_data = "replogle_rpe1_essential")), dep_jobs = list(config_job)),
          pretrained_rpe1_highDim = transfer_perturbation_prediction(c(default_params, list(ridge_penalty = 0.1, pca_dim = 40, reference_data = "replogle_rpe1_essential")), dep_jobs = list(config_job))
        ))
      }
      pert_jobs
    })
    names(inner_jobs) <- seeds
    purrr::list_flatten(inner_jobs, name_spec = "{outer}-{inner}")
  })
  names(jobs) <- datasets
  jobs <- purrr::list_flatten(jobs, name_spec = "{outer}-{inner}")
  
  params <- list(
    job_ids = map_chr(jobs, "result_id"),
    names = names(jobs)
  )
  collect_perturbation_predictions(params, dep_jobs = jobs)
}

rmap <- function(.x, .f, ..., .progress = FALSE){
  map(seq_len(nrow(.x)), \(idx){
    .f(.x[idx,], ...)
  }, .progress = .progress)
}

rmap2 <- function(.x, .y, .f, ..., .progress = FALSE){
  map2(seq_len(nrow(.x)), .y, \(idx, ..y){
    .f(.x[idx,], ..y, ...)
  }, .progress = .progress)
}

make_perturbation_sweep <- function(datasets = c('adamson', 'replogle_rpe1_essential'),
                                    ridge_penalty = 10^seq(-2, 1, length.out = 5),
                                    pca_dim = c(2, 5, 10, 20, 50, 100, 200, 500),
                                    seeds = 1:2){
  config_jobs <- tidyr::expand_grid(datasets, seeds) %>%
     mutate(config_job = rmap(across(c(dataset = datasets, seed=seeds)), prepare_perturbation_data))
  ref_jobs <- config_jobs %>%
    mutate(test_train_config_id = map_chr(config_job, "result_id")) %>%
    mutate(ground_truth = rmap2(across(c(dataset_name = datasets, test_train_config_id, seed = seeds)), config_job,
                          \(x, y) ground_truth_combinatorial_prediction(x, list(y))))

  jobs <- tidyr::expand_grid(datasets, ridge_penalty, pca_dim, seeds) %>%
    left_join(config_jobs, by = c("datasets", "seeds")) %>%
    mutate(test_train_config_id = map_chr(config_job, "result_id")) %>%
    mutate(linear = rmap2(across(c(dataset_name = datasets, test_train_config_id, ridge_penalty, pca_dim, seed = seeds)), config_job,
                         \(x, y) linear_perturbation_prediction(x, list(y))),
           pretrained = rmap2(cbind(reference_data = "replogle_k562_essential", across(c(dataset_name = datasets, test_train_config_id, ridge_penalty, pca_dim, seed = seeds))), config_job,
                          \(x, y) transfer_perturbation_prediction(x, list(y)))) %>%
    pivot_longer(c(linear, pretrained), names_to = "job_name", values_to = "job")

  sweep <- collect_perturbation_predictions(list(job_ids = map_chr(jobs$job, "result_id"), names = jobs$job_name), dep_jobs = jobs$job)
  ref <- collect_perturbation_predictions(list(job_ids = map_chr(ref_jobs$ground_truth, "result_id"), names = paste0("ground_truth-", ref_jobs$seeds, "-", ref_jobs$datasets)), dep_jobs = ref_jobs$ground_truth)
  list(sweep = sweep, ref = ref)
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

pert_jobs <- make_perturbation_jobs(seeds = 1:2)
write_rds(pert_jobs, "tmp/combinatorial_perturbation_jobs.RDS")
stat <- map_chr(pert_jobs$dependencies, job_status); table(stat)
run_job(pert_jobs, priority = "normal")
file.copy(file.path(result_file_path(pert_jobs), "predictions.RDS"), to = "output/perturbation_results_predictions.RDS")
file.copy(file.path(result_file_path(pert_jobs), "parameters.RDS"), to = "output/perturbation_results_parameters.RDS")


pert_sweep <- make_perturbation_sweep()
walk(pert_sweep, run_job, "low")


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
# Get parameter sweep output
sweep <- bind_rows(readRDS(file.path(result_file_path(pert_sweep$sweep), "predictions.RDS"))) %>%
  inner_join(readRDS(file.path(result_file_path(pert_sweep$sweep), "parameters.RDS")) %>%
               map(\(p) tibble(id = p$id, name = p$name, parameters = as_tibble(p$parameters), 
                               train = names(p$test_train_labels), perturbation = p$test_train_labels)) %>%
               bind_rows() %>%
               unnest(perturbation) %>%
               unpack(parameters),
             by = c("id", "perturbation", "name"))
ref <- bind_rows(readRDS(file.path(result_file_path(pert_sweep$ref), "predictions.RDS"))) %>%
  inner_join(readRDS(file.path(result_file_path(pert_sweep$ref), "parameters.RDS")) %>%
               map(\(p) tibble(id = p$id, name = p$name, parameters = as_tibble(p$parameters), 
                               train = names(p$test_train_labels), perturbation = p$test_train_labels)) %>%
               bind_rows() %>%
               unnest(perturbation) %>%
               unpack(parameters),
             by = c("id", "perturbation", "name")) %>%
  dplyr::select(dataset_name, test_train_config_id, seed, perturbation, obs = prediction)

sweep %>%
  tidylog::inner_join(ref, by = c("dataset_name", "test_train_config_id", "seed", "perturbation")) %>%
  tidylog::inner_join(ref %>% filter(perturbation == "ctrl") %>% dplyr::rename(baseline = obs) %>% dplyr::select(-perturbation),
                      by = c("dataset_name", "test_train_config_id", "seed")) %>%
  mutate(cor = map2_dbl(prediction, obs, cor),
         dist = map2_dbl(prediction, obs, \(x, y) sqrt(sum((x - y)^2))),
         delta_cor = pmap_dbl(list(prediction, obs, baseline), \(x, y, b) cor(x-b, y-b))) %>%
  dplyr::select(name, dataset_name, ridge_penalty, pca_dim, seed, perturbation, training, cor, dist, delta_cor) %>%
  write_tsv("output/perturbation_prediction_parameters_sweep.tsv.gz")




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


