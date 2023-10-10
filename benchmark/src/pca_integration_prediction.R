library(SingleCellExperiment)
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Run PCA integration")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--dataset_config", type = "character", nargs = 1, help = "The name of a dataset in datasets.yaml") 
pa <- argparser::add_argument(pa, "--dataset_config_override", type = "character", default = "", nargs = Inf, help = "Override settings from datasets.yaml") 

pa <- argparser::add_argument(pa, "--pca_dim", type = "integer", default = 30, nargs = 1, help = "The number of PCA dimensions")
pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The number of PCA dimensions")

pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --data_id 8f872edd6a7eb-827531ae33574
#                             --dataset_config kang
#                             --pca_dim 30
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

config <- get_data_config(pa$dataset_config, pa$dataset_config_override)

out_dir <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

sce <- zellkonverter::readH5AD(file.path(pa$working_dir, "results", pa$data_id, "train.h5ad"))
holdout <- zellkonverter::readH5AD(file.path(pa$working_dir, "results", pa$data_id, "holdout.h5ad"))

mat <- assay(sce, config$assay_continuous)
mat_holdout <- assay(sce, config$assay_continuous)
row_means <- rowMeans(mat)
pca <- irlba::prcomp_irlba(t(mat - row_means), n = pa$pca_dim, center = FALSE)


train_embedding <- t(pca$x)
holdout_embedding <- t(pca$rotation) %*% (mat_holdout - row_means)

predict_results <- replicate(length(config$contrast), NULL)
predict_holdout_results <- replicate(length(config$contrast), NULL)
names(predict_results) <- config$contrast
names(predict_holdout_results) <- config$contrast
for(lvl in config$contrast){
  predict_results[[lvl]] <- pca$rotation %*% train_embedding
  predict_holdout_results[[lvl]] <- pca$rotation %*% holdout_embedding
}

# Save everything
# qs::qsave(list(embedding = embedding), output_file)
# Save everything
dir.create(out_dir)
write.table(train_embedding, file = file.path(out_dir, glue::glue("train-embedding.tsv")), row.names = FALSE, col.names = FALSE)
write.table(holdout_embedding, file = file.path(out_dir, glue::glue("holdout-embedding.tsv")), row.names = FALSE, col.names = FALSE)

for(lvl in config$contrast){
  write.table(predict_results[[lvl]], file = file.path(out_dir, glue::glue("train-prediction_{lvl}.tsv")), row.names = FALSE, col.names = FALSE)
  write.table(predict_holdout_results[[lvl]], file = file.path(out_dir, glue::glue("holdout-prediction_{lvl}.tsv")), row.names = FALSE, col.names = FALSE)
}


#### Session Info
sessionInfo()


