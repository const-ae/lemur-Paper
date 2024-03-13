library(SingleCellExperiment)
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Run pseudobulk differential expression per cluster")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--dataset_config", type = "character", nargs = 1, help = "The name of a dataset in datasets.yaml") 
pa <- argparser::add_argument(pa, "--dataset_config_override", type = "character", default = "", nargs = Inf, help = "Override settings from datasets.yaml") 

pa <- argparser::add_argument(pa, "--pca_dim", type = "integer", default = 30, nargs = 1, help = "The number of PCA dimensions")
pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The number of PCA dimensions")

pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --data_id 35f39505c59e9-c5f2f50b80560
#                             --dataset_config kang
#                             --pca_dim 30
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

config <- get_data_config(pa$dataset_config, pa$dataset_config_override)

out_dir <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

sce <- zellkonverter::readH5AD(file.path(pa$working_dir, "results", pa$data_id, "train.h5ad"))

mat <- assay(sce, config$assay_continuous)

red_mat <- if(pa$pca_dim < min(nrow(sce), ncol(sce))){
  t(irlba::prcomp_irlba(t(mat), n = pa$pca_dim)$x)
}else{
  mat
}

embedding <- t(harmony::RunHarmony(red_mat, meta_data = colData(sce), 
                            vars_use = c(config$main_covariate, unlist(config$batch_covariates))))

# Save everything
tmp_out_dir <- paste0(out_dir, "-tmp")
dir.create(tmp_out_dir)
write.table(embedding, file = file.path(tmp_out_dir, glue::glue("train-embedding.tsv")), row.names = FALSE, col.names = FALSE, sep = "\t")
file.rename(tmp_out_dir, out_dir)

#### Session Info
sessionInfo()


