setwd("/g/huber/users/ahlmanne/projects/lemur-Paper-benchmark")
library(tidyverse)
library(glue)
library(MyWorkflowManager)
library(SingleCellExperiment)
source("submission/wrap_scripts.R")
init("/scratch/ahlmanne/lemur_benchmark")





#---------
# Make DE datasets
de_datasets <- read_tsv("data/de_data_spec.tsv")
clustering <- "kmeans"
res <- de_datasets %>%
  filter(!str_detect(vals, "\\s")) %>%
  mutate(job = map2(data, vals, \(name, cond){
    dj <- prepare_de_data(list(dataset = name, condition = cond, n_hvgs = 8000, clustering = clustering,
                               randomize = "cells", n_de_genes = 200, lfc_mean = c(0.5, 1, 2, 4), cut_at = c(2, 3, 10, 20)))
    # dj <- prepare_de_data(list(dataset = name, condition = cond, n_hvgs = 5000, clustering = clustering,
    #                            randomize = "cells", n_de_genes = 100, lfc_mean = c(1, 2, 4), cut_at = c(3, 10, 20)))
    tibble(lemur_glmGamPoi = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "glmGamPoi"), dep_job = dj)),
           lemur_edgeR = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR"), dep_job = dj)),
           lemur_limma = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "limma"), dep_job = dj)),
           lemur_edgeR_nemb_8 = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", n_embedding = 8), dep_job = dj)),
           lemur_edgeR_nemb_80 = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", n_embedding = 80), dep_job = dj)),
           lemur_edgeR_testfrac_0.8 = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", test_fraction = 0.8), dep_job = dj)),
           lemur_edgeR_skip_al = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", skip_alignment = "TRUE"), dep_job = dj)),
           lemur_edgeR_skip_multCondPCA = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", skip_multi_cond_pca = "TRUE"), dep_job = dj)),
           lemur_edgeR_dir_contrast = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", test_direction_method = "contrast"), dep_job = dj)),
           lemur_edgeR_sel_contrast = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", test_selection_method = "contrast"), dep_job = dj)),
           lemur_edgeR_sf_method_ratio = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", test_size_factor_method = "ratio"), dep_job = dj)),
           global_glmGamPoi = list(global_de(list(data_id = dj$result_id, dataset_config = name, test_method = "glmGamPoi"), dep_job = dj)),
           global_edgeR = list(global_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR"), dep_job = dj)),
           celltype_glmGamPoi = list(celltype_de(list(data_id = dj$result_id, dataset_config = name, test_method = "glmGamPoi"), dep_job = dj)),
           celltype_edgeR = list(celltype_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR"), dep_job = dj)),
           cluster_glmGamPoi = list(cluster_de(list(data_id = dj$result_id, dataset_config = name, test_method = "glmGamPoi"), dep_job = dj)),
           cluster_edgeR = list(cluster_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR"), dep_job = dj)),
           miloDE_edgeR = list(milo_de(list(data_id = dj$result_id, dataset_config = name), dep_job = dj)),
    )
  })) %>%
  unnest(job, names_sep = "-") %>%
  pivot_longer(starts_with("job"), names_sep = "-", names_to = c(".value", "method"))


res %>% mutate(status = map_chr(job, job_status)) %>% print(n = 20)
res %>% mutate(status = map_chr(job, job_status)) %>% filter(status != "done")
# walk(res$job, run_job, priority = "normal")

#---------
# Collect de results 
eval_job <- res %>%
  mutate(de_results_id = map_chr(job, "result_id"), data_results_id = map_chr(job, \(j) j$dependencies[[1]]$result_id)) %>%
  filter(!str_detect(vals, "\\s")) %>%
  {
    data <- .
    evaluate_de_results(params = c(as.list(dplyr::select(data, -job)), list(small_de_n = 400)), dep_jobs = data$job)
  }
qs::qsave(eval_job, glue("tmp/evaluation_job-{clustering}-spec.qs"))
job_status(eval_job)
run_job(eval_job, "normal")

fs::file_size(file.path(result_file_path(eval_job), "de_results.RDS"))
de_res <- readRDS(file.path(result_file_path(eval_job), "de_results.RDS"))
power_res <- readRDS(file.path(result_file_path(eval_job), "fdr_tpr_results.RDS"))

power_res %>%
  unnest(results) %>%
  dplyr::select(-strat_results) %>%
  write_tsv(glue::glue("output/differential_expression_fdr_power-{clustering}.tsv.gz"))

gene_info <- de_res %>%
  distinct(data, vals, genes) %>%
  unnest(genes)


power_per_desize <- de_res %>%
  filter(str_detect(method, "edgeR")) %>%
  mutate(de_truth = map2(de, ground_truth, \(de, gt){
    d <- if("cell_type" %in% colnames(de)){
      inner_join(de, gt, by = c("name" = "gene", "cell_type" = "group"))
    }else{
      inner_join(de, gt, by = c("name" = "gene"))
    }
    if("n_cells" %in% colnames(d)){
      d <- dplyr::rename(d, nei_n_cells = n_cells)
    }
    d
  })) %>%
  dplyr::select(-c(de, ground_truth)) %>%
  unnest(de_truth) %>%
  tidylog::left_join(gene_info, by = c("data", "vals", "name")) %>%
  filter(de_size != 0) %>% 
  filter(! is.na(is_modified) & is_modified) %>%
  mutate(data_val = paste0(data, "-", vals)) %>%
  mutate(data_val = fct_reorder(data_val, n_cells)) %>%
  mutate(adj_pval = ifelse(is.na(adj_pval), 1, adj_pval)) %>%
  summarize(signif = any(adj_pval < 0.1), de_size = mean(de_size), .by = c(method, data_val, name)) 

power_per_desize %>%
  write_tsv(glue::glue("output/differential_expression_fdr_power-{clustering}-stratified.tsv.gz"))




#------------

cell_names <- res %>% 
  mutate(data_job = map(job, \(x) x$dependencies[[1]])) %>%
  distinct(data, vals, data_job) %>%
  mutate(cell_names = map(data_job, \(dj){
    sce <- qs::qread(result_file_path(dj))
    colnames(sce)
  }))


lemur_recall_prec <- de_res %>%
  filter(method == "lemur_edgeR") %>%
  mutate(de_truth = map2(de, ground_truth, \(de, gt){
    d <- if("cell_type" %in% colnames(de)){
      inner_join(de, gt, by = c("name" = "gene", "cell_type" = "group"))
    }else{
      inner_join(de, gt, by = c("name" = "gene"))
    }
    if("n_cells" %in% colnames(d)){
      d <- dplyr::rename(d, nei_n_cells = n_cells)
    }
    d
  })) %>%
  dplyr::select(-c(de, ground_truth)) %>%
  unnest(de_truth) %>%
  tidylog::left_join(gene_info, by = c("data", "vals", "name")) %>%
  filter(de_size != 0) %>%
  dplyr::select(data, vals, name, is_de_cell, neighborhood, de_n_cells = n_cells, nei_n_cells, adj_pval)

lemur_nei_overlap <- lemur_recall_prec %>%
  # filter(adj_pval < 0.1) %>%
  left_join(cell_names, by = c("data", "vals")) %>%
  mutate(overlap = pmap(list(is_de_cell, cell_names, neighborhood), \(is_de, names, nei){
    true_changed <- names[is_de]
    tibble(TP = length(intersect(true_changed, nei)), FP = length(nei) - TP,
           FN = length(true_changed) - TP, de_n_cells = length(true_changed), nei_size = length(nei))
  })) %>%
  dplyr::select(data, vals, name, adj_pval, overlap) %>%
  unnest(overlap)

lemur_nei_overlap %>%
  write_tsv("output/differential_expression-kang-recall_precision.tsv.gz")

lemur_nei_overlap %>%
  mutate(recall = TP / de_n_cells, 
         precision = TP / nei_size) %>%
  ggplot(aes(x = recall, y = precision)) +
    geom_point(aes(color = adj_pval < 0.1))

#------------
kang_lemur_res <- res %>%
  filter(data == "kang" & vals == "ctrl" & method == "lemur_edgeR") %>%
  mutate(data_job = map(job, \(d) d$dependencies[[1]]))

wd <- MyWorkflowManager:::.OUTPUT_FOLDER()
sce <- qs::qread(file.path(wd, "results", kang_lemur_res$data_job[[1]]$result_id))
pa <- list(n_embedding = 30, test_fraction = 0.5, test_method = "edgeR")
library(lemur)
set.seed(1)
fit <- lemur(sce, design = ~ fake_condition, n_embedding = pa$n_embedding, 
             use_assay = "logcounts", test_fraction = pa$test_fraction)
fit <- align_harmony(fit)
fit <- test_de(fit, contrast = cond(fake_condition = "fake_trt") - cond(fake_condition = "fake_ctrl"))
nei <- find_de_neighborhoods(fit, group_by = vars(fake_condition, sample), test_method = pa$test_method, count_assay_name = "counts")

nei_alt <- readRDS(file.path(result_file_path(kang_lemur_res$job[[1]]), "de_results.RDS"))
testthat::expect_equal(nei, nei_alt)

saveRDS(fit[rowData(fit)$is_simulated,], "output/differential-expression-kang_lemur_fit.RDS")
saveRDS(nei[rowData(fit)$is_simulated,], "output/differential-expression-kang_lemur_fit-neighborhood.RDS")

