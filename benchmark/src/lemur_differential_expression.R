library(SingleCellExperiment)
library(lemur)
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Run lemur")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--dataset_config", type = "character", nargs = 1, help = "The name of a dataset in datasets.yaml") 
pa <- argparser::add_argument(pa, "--dataset_config_override", type = "character", default = "", nargs = Inf, help = "Override settings from datasets.yaml") 

pa <- argparser::add_argument(pa, "--n_embedding", type = "integer", default = 30, nargs = 1, help = "The number of PCA dimensions")
pa <- argparser::add_argument(pa, "--test_fraction", type = "numeric", default = 0.5, nargs = 1, help = "The fraction of cells used for testing")
pa <- argparser::add_argument(pa, "--skip_multi_cond_pca", type = "character", default = "FALSE", help = "Set to true, to fix design to ~ 1 in first lemur fit")
pa <- argparser::add_argument(pa, "--skip_alignment", type = "character", default = "FALSE", help = "Set to true, to skip harmony")
pa <- argparser::add_argument(pa, "--test_method", type = "character", default = "glmGamPoi", help = "Select the test method")
pa <- argparser::add_argument(pa, "--test_direction_method", type = "character", default = "random", help = "Select the test method")
pa <- argparser::add_argument(pa, "--test_selection_method", type = "character", default = "zscore", help = "Select the test method")
pa <- argparser::add_argument(pa, "--test_size_factor_method", type = "character", default = "normed_sum", help = "Select the test method")


pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The seed")


pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --data_id a7d141f51613f-cabcc10d4ce8e
#                             --dataset_config angelidis
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

config <- get_data_config(pa$dataset_config, pa$dataset_config_override)
skip_multi_cond_pca <- as.logical(yaml::read_yaml(text = pa$skip_multi_cond_pca))
skip_alignment <- as.logical(yaml::read_yaml(text = pa$skip_alignment))
out_dir <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

sce <- qs::qread(file.path(pa$working_dir, "results", pa$data_id))

batch_covs <- paste0(vapply(config$batch_covariates, \(x) paste0(" + ", x), FUN.VALUE = character(1L)), collapse = "", recycle0 = TRUE)
full_formula_string <-  paste0("~ fake_condition ", batch_covs)
formula_string <- if(skip_multi_cond_pca){ 
  "~ 1"
}else{ 
  full_formula_string
}
lemur_fit_time <- system.time({
  fit <- lemur(sce, design = as.formula(formula_string), n_embedding = pa$n_embedding, 
                      use_assay = config$assay_continuous, test_fraction = pa$test_fraction)
})

lemur_align_time <- system.time({
  if(! skip_alignment){
    fit <- align_harmony(fit, design = as.formula(full_formula_string))
  }
})

if(skip_multi_cond_pca){
  lemur_test_time <- system.time({
    pred1 <- predict(fit, newdesign = 1, alignment_design_matrix =  lemur:::parse_contrast(cond(fake_condition = "fake_trt"), fit$alignment_design, simplify = TRUE))
    pred2 <- predict(fit, newdesign = 1, alignment_design_matrix =  lemur:::parse_contrast(cond(fake_condition = "fake_ctrl"), fit$alignment_design, simplify = TRUE))
  })
  lemur_nei_time <- system.time({
    nei <- find_de_neighborhoods(fit, group_by = vars(fake_condition, sample), test_method = pa$test_method,
                                 design = fit$alignment_design, contrast = cond(fake_condition = "fake_trt") - cond(fake_condition = "fake_ctrl"), 
                                 de_mat = pred1 - pred2,
                                 selection_procedure = pa$test_selection_method,
                                 directions = pa$test_direction_method, size_factor_method = pa$test_size_factor_method,
                                 count_assay_name = config$assay_counts)                                                      
  })
  
}else{
  lemur_test_time <- system.time({
    fit <- test_de(fit, contrast = cond(fake_condition = "fake_trt") - cond(fake_condition = "fake_ctrl"))
  })
  lemur_nei_time <- system.time({
    nei <- find_de_neighborhoods(fit, group_by = vars(fake_condition, sample), test_method = pa$test_method,
                                 selection_procedure = pa$test_selection_method,
                                 directions = pa$test_direction_method, size_factor_method = pa$test_size_factor_method,
                                 count_assay_name = config$assay_counts)
  })
}

print(table(nei$adj_pval < 0.1, fit$rowData$is_simulated))

# Save everything
dir.create(out_dir)
saveRDS(nei, file.path(out_dir, "de_results.RDS"))
saveRDS(list(lemur_fit_time = lemur_fit_time, lemur_align_time = lemur_align_time,
               lemur_test_time = lemur_test_time, lemur_nei_time = lemur_nei_time), file.path(out_dir, "duration.RDS"))

#### Session Info
sessionInfo()


