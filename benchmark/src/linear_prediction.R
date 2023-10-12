library(SingleCellExperiment)
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Run pseudobulk differential expression per cluster")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--dataset_config", type = "character", nargs = 1, help = "The name of a dataset in datasets.yaml") 
pa <- argparser::add_argument(pa, "--dataset_config_override", type = "character", default = "", nargs = Inf, help = "Override settings from datasets.yaml") 

pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --data_id 8f872edd6a7eb-827531ae33574
#                             --dataset_config kang
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)

config <- get_data_config(pa$dataset_config, pa$dataset_config_override)

out_dir <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

sce <- zellkonverter::readH5AD(file.path(pa$working_dir, "results", pa$data_id, "train.h5ad"))
holdout <- zellkonverter::readH5AD(file.path(pa$working_dir, "results", pa$data_id, "holdout.h5ad"))

mat <- as.matrix(assay(sce, config$assay_continuous))
batch_covs <- paste0(vapply(config$batch_covariates, \(x) paste0(" + ", x), FUN.VALUE = character(1L)), collapse = "", recycle0 = TRUE)
formula_string <-  paste0("t(mat) ~ ", config$main_covariate, " ", batch_covs)
df <- colData(sce)
df[[config$main_covariate]] <- as.factor(df[[config$main_covariate]])
df_holdout <- colData(holdout)
df_holdout[[config$main_covariate]] <- factor(df_holdout[[config$main_covariate]], 
                                              levels = levels(df[[config$main_covariate]]))
  
fit <- lm(as.formula(formula_string), data = df)

predict_results <- replicate(length(config$contrast), NULL)
predict_holdout_results <- replicate(length(config$contrast), NULL)
names(predict_results) <- config$contrast
names(predict_holdout_results) <- config$contrast
for(lvl in config$contrast){
  df[[config$main_covariate]] <- lvl
  df_holdout[[config$main_covariate]] <- lvl
  predict_results[[lvl]] <- t(predict(fit, df))
  predict_holdout_results[[lvl]] <- t(predict(fit, df_holdout))
}


# # Save everything
# qs::qsave(list(predictions = predict_results), output_file)
# Save everything
dir.create(out_dir)
for(lvl in config$contrast){
  write.table(predict_results[[lvl]], file = file.path(out_dir, glue::glue("train-prediction_{lvl}.tsv")), row.names = FALSE, col.names = FALSE, sep = "\t")
  write.table(predict_holdout_results[[lvl]], file = file.path(out_dir, glue::glue("holdout-prediction_{lvl}.tsv")), row.names = FALSE, col.names = FALSE, sep = "\t")
}

#### Session Info
sessionInfo()


