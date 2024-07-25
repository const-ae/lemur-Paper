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
  mutate(seed = list(1:10)) %>%
  unnest(seed) %>%
  filter(!str_detect(vals, "\\s")) %>%
  mutate(job = pmap(list(data, vals, seed), \(name, cond, se){
    dj <- prepare_de_data(list(dataset = name, condition = cond, n_hvgs = 8000, clustering = clustering, seed  = se,
                               randomize = "cells", n_de_genes = 200, lfc_mean = c(0.5, 1, 2, 4), cut_at = c(2, 3, 10, 20)))
    # dj <- prepare_de_data(list(dataset = name, condition = cond, n_hvgs = 5000, clustering = clustering,
    #                            randomize = "cells", n_de_genes = 100, lfc_mean = c(1, 2, 4), cut_at = c(3, 10, 20)))
    tibble(lemur_glmGamPoi = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "glmGamPoi"), dep_job = dj)),
           lemur_edgeR = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR"), dep_job = dj)),
           lemur_limma = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "limma"), dep_job = dj)),
           lemur_edgeR_notesttraining = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", split_test_training = "FALSE"), dep_job = dj)),
           lemur_edgeR_count_split = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", do_count_splitting = "TRUE"), dep_job = dj)),
           lemur_edgeR_nemb_8 = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", n_embedding = 8), dep_job = dj)),
           lemur_edgeR_nemb_80 = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", n_embedding = 80), dep_job = dj)),
           lemur_edgeR_testfrac_0.8 = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", test_fraction = 0.2), dep_job = dj)),
           lemur_edgeR_testfrac_0.2 = list(lemur_de(list(data_id = dj$result_id, dataset_config = name, test_method = "edgeR", test_fraction = 0.8), dep_job = dj)),
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
# walk(res$job, run_job, priority = "low")
stat <- map_chr(res$job, job_status); table(stat)

#---------
# Collect de results 
eval_jobs <- res %>%
  filter(!str_detect(vals, "\\s")) %>%
  mutate(split_key = rep_len(1:30, length.out = n())) %>%
  dplyr::rename(data_seeds = seed) %>%
  mutate(de_results_id = map_chr(job, "result_id"), data_results_id = map_chr(job, \(j) j$dependencies[[1]]$result_id)) %>%
  group_by(split_key) %>%
  group_map(\(data, key){
    evaluate_de_results(params = c(as.list(dplyr::select(data, -job)), list(small_de_n = 400)), dep_jobs = data$job, 
                        memory = "80GB", duration = "02:00:00")
  })
qs::qsave(eval_jobs, glue("tmp/evaluation_job-{clustering}-spec.qs"))
table(map_chr(eval_jobs, job_status))
walk(eval_jobs, \(j) run_job(j, "normal"))

power_res <- bind_rows(map(eval_jobs, \(j){
  readRDS(file.path(result_file_path(j), "fdr_tpr_results.RDS"))  
}))

power_res %>%
  unnest(results) %>%
  dplyr::select(-strat_results) %>%
  write_tsv(glue::glue("output/differential_expression_fdr_power-{clustering}.tsv.gz"))



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

