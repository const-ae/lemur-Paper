library(SingleCellExperiment)
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Run pseudobulk differential expression per cluster")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--dataset_config", type = "character", nargs = 1, help = "The name of a dataset in datasets.yaml") 
pa <- argparser::add_argument(pa, "--dataset_config_override", type = "character", default = "", nargs = Inf, help = "Override settings from datasets.yaml") 

pa <- argparser::add_argument(pa, "--latent_dim", type = "integer", default = 30, nargs = 1, help = "The number of PCA dimensions")
pa <- argparser::add_argument(pa, "--skip_multi_cond_pca", type = "flag", default = FALSE, help = "Skip the multi-condition PCA step")
pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The number of PCA dimensions")

pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --data_id 35f39505c59e9-c5f2f50b80560
#                             --dataset_config kang
#                             --latent_dim 30
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

config <- get_data_config(pa$dataset_config, pa$dataset_config_override)

output_file <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

sce <- qs::qread(file.path(pa$working_dir, "results", pa$data_id))$train
sce$split <- sample(c("training", "valid", "ood"), size = ncol(sce), prob = c(0.8, 0.1, 0.1), replace = TRUE)

adata_file <- tempfile()
embedding_out_file <- tempfile()

hvg <- order(-rowVars(assay(sce, "logcounts")))
sce_subset <- sce[hvg[1:1000], sample.int(ncol(sce), size = 5000)]
zellkonverter::writeH5AD(sce_subset, adata_file, X_name = "logcounts")
sce2 <- zellkonverter::readH5AD(adata_file)

tmp <- assays(sce2)[[1]]
assays(sce2)[[1]] <- assays(sce2)[[2]]
assays(sce2)[[2]] <- tmp
colData(sce2)
colData(sce_subset)
testthat::expect_equal(sce_subset, sce2)

# Call CPA
launch_command <- glue::glue(r"(
src/python_scripts/python_script_wrapper.sh cpa_env \
  src/python_scripts/run_cpa.py \
    --adata_file {adata_file} \
    --embedding_out_file {embedding_out_file} \
    --perturbation_key {config$main_covariate} \
    --split_key split \
    --control_group {config$contrast[1]} \
    --n_latent 15 \
    --max_epochs 100
)")
message(launch_command)
system(launch_command)

latent_after <- t(as.matrix(read.delim(embedding_out_file, header = FALSE, sep = "\t")))
umap_after <- uwot::umap(t(latent_after))
as_tibble(colData(sce_subset)) %>%
  mutate(umap = umap_after) %>%
  ggplot(aes(x = umap[,1], y = umap[,2])) +
    geom_point(aes(color = label))

# Save everything
qs::qsave(list(embedding = embedding), output_file)


#### Session Info
sessionInfo()


