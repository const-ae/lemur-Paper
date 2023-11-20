library(SingleCellExperiment)
source("src/utils/data_helper.R")
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Run pseudobulk differential expression per cluster")
pa <- argparser::add_argument(pa, "--dataset", type = "character", nargs = 1, help = "The name of a dataset") 
pa <- argparser::add_argument(pa, "--variation", type = "character", nargs = 1, help = "The data variation that is generated") 
pa <- argparser::add_argument(pa, "--seed", type = "integer", nargs = 1, default = 1, help = "Seed to tame randomness") 
pa <- argparser::add_argument(pa, "--hold_out_percentage", type = "", nargs = "numeric", default = 0.2, help = "The holdout percentage if variation is 'random_holdout'") 
pa <- argparser::add_argument(pa, "--n_hvgs", type = "", nargs = "numeric", default = 500, help = "The number of genes selected as highly variable") 
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --dataset kang
#                             --variation identity
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

if(! pa$dataset %in% names(data_loader)){
  stop("No data loader defined for ", pa$dataset, ".\nThe options are ", toString(names(data_loader)))
}

valid_variations <- c("random_holdout", "random_holdout_hvg", "cluster_holdout", "cluster_holdout_hvg", "identity", "small")
var <- pa$variation
sel_hvg <- var %in% c("small", "random_holdout_hvg", "cluster_holdout_hvg")
if(! var %in% valid_variations){
  stop("Variation ", var,  " is not a recognized dataset variation")
}

output_file <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

config <- get_data_config(pa$dataset)
sce <- data_loader[[pa$dataset]]()

if(sel_hvg){
  hvg <- order(-rowVars(assay(sce, config$assay_continuous)))
  n_hvgs <- min(nrow(sce), pa$n_hvgs)
  sce <- sce[hvg[seq_len(n_hvgs)],]
}

if(var == "identity"){
  train <- sce
  holdout <- NULL
}else if(var %in% c("random_holdout", "random_holdout_hvg")){
  holdout_lgl <- runif(seq_len(ncol(sce))) < pa$hold_out_percentage
  train <- sce[,!holdout_lgl]
  holdout <- sce[,holdout_lgl]
}else if(var %in% c("cluster_holdout", "cluster_holdout_hvg")){
  data <- assay(sce, config$assay_continuous)
  is_gr1 <- colData(sce)[[config$main_covariate]] == config$contrast[1]
  is_gr2 <- colData(sce)[[config$main_covariate]] == config$contrast[2]
  pca_ref_gr1 <- irlba::prcomp_irlba(t(data[, is_gr1]), n = 30)  
  pca_ref_gr2 <- irlba::prcomp_irlba(t(data[, is_gr2]), n = 30)  
  clustering_ref_gr1 <- bluster::clusterRows(pca_ref_gr1$x, bluster::KNNGraphParam())
  clustering_ref_gr2 <- bluster::clusterRows(pca_ref_gr2$x, bluster::KNNGraphParam())
  holdout_lgl <- rep(FALSE, ncol(sce))
  holdout_lgl[is_gr1] <- clustering_ref_gr1 == "1"
  holdout_lgl[is_gr2] <- clustering_ref_gr2 == "1"
  
  train <- sce[,!holdout_lgl]
  holdout <- sce[,holdout_lgl]  
}else if(var == "small"){
  n_cells <- min(ncol(sce), 5000)
  train <- sce[, sample.int(ncol(sce), size = n_cells)]
  holdout <- sce[, sample.int(ncol(sce), size = n_cells)]
}

# Save everything
if(dir.exists(output_file)){
  unlink(output_file, recursive = TRUE)
}
dir.create(output_file)
zellkonverter::writeH5AD(train, file.path(output_file, "train.h5ad"), X_name = config$assay_continuous)
if(! is.null(holdout)){
  zellkonverter::writeH5AD(holdout, file.path(output_file, "holdout.h5ad"), X_name = config$assay_continuous)
}

# Session Info
sessionInfo()