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
#                             --data_id 35f39505c59e9-c5f2f50b8056
#                             --dataset_config kang
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

config <- get_data_config(pa$dataset_config, pa$dataset_config_override)

output_file <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

sce <- qs::qread(file.path(pa$working_dir, "results", pa$data_id))$train

mat <- as.matrix(assay(sce, config$assay_continuous))
batch_covs <- paste0(vapply(config$batch_covariates, \(x) paste0(" + ", x), FUN.VALUE = character(1L)), collapse = "", recycle0 = TRUE)
formula_string <-  paste0("t(mat) ~ ", config$main_covariate, " ", batch_covs)
df <- colData(sce)
df[[config$main_covariate]] <- as.factor(df[[config$main_covariate]])

fit <- lm(as.formula(formula_string), data = df)

predict_results <- replicate(length(config$contrast), NULL)
names(predict_results) <- config$contrast
for(lvl in config$contrast){
  df[[config$main_covariate]] <- lvl
  predict_results[[lvl]] <- t(predict(fit, df))
}


# Save everything
qs::qsave(list(predictions = predict_results), output_file)


#### Session Info
sessionInfo()


