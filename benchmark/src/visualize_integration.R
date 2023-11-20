library(SingleCellExperiment)
library(tidyverse)
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Prepare visualization of the integration results")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--dataset_config", type = "character", nargs = 1, help = "The name of a dataset in datasets.yaml") 
pa <- argparser::add_argument(pa, "--dataset_config_override", type = "character", default = "", nargs = Inf, help = "Override settings from datasets.yaml") 
pa <- argparser::add_argument(pa, "--integration_id", type = "character", help = "The result id of the integration job")

pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The seed")

pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --data_id 8f872edd6a7eb-3c3d82e59c413
#                             --dataset_config canogamez
#                             --integration_id 8e4417d74d625-706dcf1809e03
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

config <- get_data_config(pa$dataset_config, pa$dataset_config_override)
outfile <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

sce <- zellkonverter::readH5AD(file.path(pa$working_dir, "results", pa$data_id, "train.h5ad"))
holdout <- zellkonverter::readH5AD(file.path(pa$working_dir, "results", pa$data_id, "holdout.h5ad"))

train_file <- file.path(pa$working_dir, "results", pa$integration_id, "train-embedding.tsv")
holdout_file <- file.path(pa$working_dir, "results", pa$integration_id, "holdout-embedding.tsv")
ref_embedding <- as.matrix(data.table::fread(train_file, sep = "\t", header = FALSE, col.names = colnames(sce)))

ref_col_data <- tibble(name = colnames(sce), 
                       covariate = colData(sce)[[config$main_covariate]],
                       sample = colData(sce)[[config$sample_covariate]],
                       celltype = colData(sce)[[config$cell_type_column]])
if(length(config$batch_covariates) > 0){
  ref_col_data <- bind_cols(ref_col_data, as_tibble(colData(sce)[config$batch_covariates]))
}


if(file.exists(holdout_file)){
  holdout_embedding <- as.matrix(data.table::fread(holdout_file, sep = "\t", header = FALSE, col.names = colnames(holdout)))
  holdout_col_data <- tibble(name = colnames(holdout), 
                         covariate = colData(holdout)[[config$main_covariate]],
                         sample = colData(holdout)[[config$sample_covariate]],
                         celltype = colData(holdout)[[config$cell_type_column]])
  if(length(config$batch_covariates) > 0){
    holdout_col_data <- bind_cols(holdout_col_data, as_tibble(colData(holdout)[config$batch_covariates]))
  }
  
  res <- rbind(ref_col_data, holdout_col_data)
  res$is_holdout <- rep(c(FALSE, TRUE), times = c(ncol(sce), ncol(holdout)))
  embedding <- cbind(ref_embedding, holdout_embedding)
}else{
  res <- ref_col_data
  embedding <- ref_embedding
}

umap <- uwot::umap(t(embedding))
tsne <- snifter::fitsne(t(embedding))
pca <- if(nrow(embedding) <= 10){
  unname(prcomp(t(embedding), rank. = 2)$x)
}else{
  unname(irlba::prcomp_irlba(t(embedding), n = 2)$x)
}

res$umap1 <- umap[,1]
res$umap2 <- umap[,2]
res$tsne1 <- tsne[,1]
res$tsne2 <- tsne[,2]
res$pca1 <- pca[,1]
res$pca2 <- pca[,2]

readr::write_tsv(res, file = outfile)


# Session Info
sessionInfo()
