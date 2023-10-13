library(SingleCellExperiment)
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Run lemur")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--dataset_config", type = "character", nargs = 1, help = "The name of a dataset in datasets.yaml") 
pa <- argparser::add_argument(pa, "--dataset_config_override", type = "character", default = "", nargs = Inf, help = "Override settings from datasets.yaml") 

pa <- argparser::add_argument(pa, "--n_embedding", type = "integer", default = 30, nargs = 1, help = "The number of PCA dimensions")
pa <- argparser::add_argument(pa, "--skip_multi_cond_pca", type = "character", default = "FALSE", help = "Set to true, to fix design to ~ 1 in first lemur fit")
pa <- argparser::add_argument(pa, "--skip_alignment", type = "character", default = "FALSE", help = "Set to true, to skip harmony")
pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The seed")


pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --data_id 8f872edd6a7eb-827531ae33574
#                             --dataset_config kang
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

config <- get_data_config(pa$dataset_config, pa$dataset_config_override)
skip_multi_cond_pca <- as.logical(yaml::read_yaml(text = pa$skip_multi_cond_pca))
skip_alignment <- as.logical(yaml::read_yaml(text = pa$skip_alignment))
out_dir <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

sce <- zellkonverter::readH5AD(file.path(pa$working_dir, "results", pa$data_id, "train.h5ad"))
holdout <- zellkonverter::readH5AD(file.path(pa$working_dir, "results", pa$data_id, "holdout.h5ad"))


batch_covs <- paste0(vapply(config$batch_covariates, \(x) paste0(" + ", x), FUN.VALUE = character(1L)), collapse = "", recycle0 = TRUE)
full_formula_string <-  paste0("~ ", config$main_covariate, " ", batch_covs)
formula_string <- if(skip_multi_cond_pca){ 
  "~ 1"
}else{ 
  full_formula_string
}
fit <- lemur::lemur(sce, design = as.formula(formula_string), n_embedding = pa$n_embedding, use_assay = config$assay_continuous,
                    test_fraction = 0)
if(! skip_alignment){
  fit <- lemur::align_harmony(fit, design = as.formula(full_formula_string))
}

train_embedding <- fit$embedding
holdout_fit <- lemur::project_on_lemur_fit(fit, holdout, return = "lemur_fit")
holdout_embedding <- holdout_fit$embedding


# Predict all conditions
predict_results <- replicate(length(config$contrast), NULL)
predict_holdout_results <- replicate(length(config$contrast), NULL)
if(skip_multi_cond_pca){
  for(idx in seq_along(config$contrast)){
    lvl <- config$contrast[idx]
    des <- lemur:::parse_contrast(cond(!!config$main_covariate := lvl), fit$alignment_design)  
    alignment_design_matrix <- matrix(des, nrow = ncol(fit), ncol = length(des))
    alignment_design_matrix2 <- matrix(des, nrow = ncol(holdout_fit), ncol = length(des))
    predict_results[[idx]] <- predict(fit, alignment_design_matrix = alignment_design_matrix)
    predict_holdout_results[[idx]] <- predict(holdout_fit, alignment_design_matrix = alignment_design_matrix2)
  }
}else{
  for(idx in seq_along(config$contrast)){
    lvl <- config$contrast[idx]
    predict_results[[idx]] <- predict(fit, newcondition = cond(!!config$main_covariate := lvl))
    predict_holdout_results[[idx]] <- predict(holdout_fit, newcondition = cond(!!config$main_covariate := lvl))
  }
}

# Save everything
dir.create(out_dir)
write.table(train_embedding, file = file.path(out_dir, glue::glue("train-embedding.tsv")), row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(holdout_embedding, file = file.path(out_dir, glue::glue("holdout-embedding.tsv")), row.names = FALSE, col.names = FALSE, sep = "\t")
for(idx in seq_along(config$contrast)){
  write.table(predict_results[[idx]], file = file.path(out_dir, glue::glue("train-prediction_{config$contrast[idx]}.tsv")), row.names = FALSE, col.names = FALSE, sep = "\t")
  write.table(predict_holdout_results[[idx]], file = file.path(out_dir, glue::glue("holdout-prediction_{config$contrast[idx]}.tsv")), row.names = FALSE, col.names = FALSE, sep = "\t")
}


#### Session Info
sessionInfo()


