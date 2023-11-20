library(SingleCellExperiment)
library(tidyverse)
source("src/utils/config_helper.R")

rbf_kernel <- function(x, y, gamma = NULL){
  if(is.null(gamma)){
    gamma <- 1/ncol(x)
  }
  exp(-gamma * rdist::cdist(x, y)^2)
}

#' Maximum mean discrepancy
#' 
#' Implemented based on cellot/losses/mmd
#' 
#' @param x,y matrix with observations in the rows and 
#'   features in the columns Note: This is the transpose of the usual data
#'   storage!
#' @param gamma the scale of the RBF kernel. If `NULL` it is `1/ncol(x)`.
#' 
mmd_distance <- function(x, y, gamma){
  xx <- rbf_kernel(x, x, gamma)
  yy <- rbf_kernel(y, y, gamma)
  xy <- rbf_kernel(x, y, gamma)
  
  mean(xx) + mean(yy) - 2 * mean(xy)
}


mmd_distance_wrapper <- function(data, ref_data, gammas = 10^seq(1, -3, length.out = 50), n_cells = c(100, 200, 500)){
  if(nrow(data) > 100){
    message("Reduce number of features to 100 for MMD.")
    all_data <- cbind(data, ref_data)
    centers <- rowMeans(all_data)
    pca <- irlba::prcomp_irlba(t(all_data), n = 50, center = TRUE)
    data <- t(pca$rotation) %*% (data - centers)
    ref_data <- t(pca$rotation) %*% (ref_data - centers)
  }
  
  
  res <- rep(NA, length(gammas) * length(n_cells))
  idx <- 1
  
  for(nc in n_cells){
    if(nc > min(ncol(data), ncol(ref_data))){
      break
    }
    sel_e <- sample.int(ncol(data), size = nc)
    sel_r <- sample.int(ncol(ref_data), size = nc)
    datat_sub <- t(data[,sel_e,drop=FALSE])
    ref_datat_sub <- t(ref_data[,sel_r,drop=FALSE])
    for(g in gammas){
      res[idx] <- mmd_distance(datat_sub, ref_datat_sub, gamma = g)
      idx <- idx + 1
    }
  }
  mean(res, na.rm = TRUE)
}


wasserstein_wrapper <- function(data, ref_data, n_cells = c(100, 200, 500)){
  res <- rep(NA, length(n_cells))
  idx <- 1
  
  for(nc in n_cells){
    if(nc > min(ncol(data), ncol(ref_data))){
      break
    }
    sel_e <- sample.int(ncol(data), size = nc)
    sel_r <- sample.int(ncol(ref_data), size = nc)
    datat_sub <- t(data[,sel_e,drop=FALSE])
    ref_datat_sub <- t(ref_data[,sel_r,drop=FALSE])
    
    tryCatch({
      res[idx] <- transport::wasserstein(transport::pp(datat_sub), transport::pp(ref_datat_sub))
    }, error = function(err){
      message("Optimal transport failed for ", nc)
    })
    idx <- idx + 1
  }
  mean(res, na.rm = TRUE)
}


mean_sd_agreement <- function(data, ref_data){
  em <- sparseMatrixStats::rowMeans2(data, useNames = FALSE)
  rm <- sparseMatrixStats::rowMeans2(ref_data, useNames = FALSE)
  
  es <- sparseMatrixStats::rowSds(data, useNames = FALSE)
  rs <- sparseMatrixStats::rowSds(ref_data, useNames = FALSE)
  
  l2_dist <- function(x, y) sqrt(sum((x - y)^2))
  
  list(l2_mean = l2_dist(rm, em), l2_sd = l2_dist(rs, es), r2_mean = cor(rm, em), r2_sd = cor(rs, es))
}

mean_sd_agreement_per_celltype <- function(data, ref_data, cell_types, ref_cell_types){
  res <- as.list(colMeans(do.call(rbind, lapply(as.character(unique(cell_types)), \(ct){
    unlist(mean_sd_agreement(data[,cell_types == ct,drop=FALSE], ref_data[,ref_cell_types == ct,drop=FALSE]))
  })), na.rm = TRUE))
  names(res) <- paste0(names(res), "_per_celltype")
  res
}

calculate_metrics <- function(data, ref_data, cell_types, ref_cell_types){
  purrr::flatten(list(
    mean_sd_agreement(data, ref_data),
    mean_sd_agreement_per_celltype(data, ref_data, cell_types, ref_cell_types),
    mmd = mmd_distance_wrapper(data, ref_data),
    wasserstein = wasserstein_wrapper(data, ref_data)
  ))
}


pa <- argparser::arg_parser("Evaluate prediction results")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--dataset_config", type = "character", nargs = 1, help = "The name of a dataset in datasets.yaml") 
pa <- argparser::add_argument(pa, "--dataset_config_override", type = "character", default = "", nargs = Inf, help = "Override settings from datasets.yaml") 
pa <- argparser::add_argument(pa, "--prediction_id", type = "character", help = "The result id of the prediction job")

pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The seed")

pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --data_id f1809110732c2-e5f290978d2e2
#                             --dataset_config kang
#                             --prediction_id 8e4417d74d625-1e9271d3f54d6
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

config <- get_data_config(pa$dataset_config, pa$dataset_config_override)
outfile <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

sce <- zellkonverter::readH5AD(file.path(pa$working_dir, "results", pa$data_id, "train.h5ad"))
holdout <- zellkonverter::readH5AD(file.path(pa$working_dir, "results", pa$data_id, "holdout.h5ad"))

train_file1 <- file.path(pa$working_dir, "results", pa$prediction_id, paste0("train-prediction_", config$contrast[1], ".tsv"))
train_file2 <- file.path(pa$working_dir, "results", pa$prediction_id, paste0("train-prediction_", config$contrast[2], ".tsv"))
holdout_file1 <- file.path(pa$working_dir, "results", pa$prediction_id, paste0("holdout-prediction_", config$contrast[1], ".tsv"))
holdout_file2 <- file.path(pa$working_dir, "results", pa$prediction_id, paste0("holdout-prediction_", config$contrast[2], ".tsv"))

ref_is_gr1 <- colData(sce)[[config$main_covariate]] == config$contrast[1]
ref_is_gr2 <- colData(sce)[[config$main_covariate]] == config$contrast[2]
ref_obs <- as.matrix(assay(sce, config$assay_continuous))

ref_obs1 <- ref_obs[,ref_is_gr1]
ref_obs2 <- ref_obs[,ref_is_gr2]
ref_celltypes <- colData(sce)[[config$cell_type_column]]


# Against training
train_pred1 <- as.matrix(data.table::fread(train_file1, sep = "\t", header = FALSE, col.names = colnames(sce)))
train_pred2 <- as.matrix(data.table::fread(train_file2, sep = "\t", header = FALSE, col.names = colnames(sce)))
## Subset pred to cross-condition (i.e., cell is ctrl and cond is trt)
train_pred1 <- train_pred1[,!ref_is_gr1]
train_pred2 <- train_pred2[,!ref_is_gr2]
res_train1 <- calculate_metrics(train_pred1, ref_obs1, ref_celltypes[!ref_is_gr1], ref_celltypes[ref_is_gr1])
res_train1$comparison <- "train1_vs_obs1"
res_train2 <- calculate_metrics(train_pred2, ref_obs2, ref_celltypes[!ref_is_gr2], ref_celltypes[ref_is_gr2])
res_train2$comparison <- "train2_vs_obs2"


res <- if(all(file.exists(c(holdout_file1, holdout_file2)))){
  pred1 <- as.matrix(data.table::fread(holdout_file1, sep = "\t", header = FALSE))
  pred2 <- as.matrix(data.table::fread(holdout_file2, sep = "\t", header = FALSE))
  is_gr1 <- colData(holdout)[[config$main_covariate]] == config$contrast[1]
  is_gr2 <- colData(holdout)[[config$main_covariate]] == config$contrast[2]
  celltypes <- colData(holdout)[[config$cell_type_column]]
  pred1 <- pred1[,!is_gr1]
  pred2 <- pred2[,!is_gr2]
  
  res_holdout1 <- calculate_metrics(pred1, ref_obs1, celltypes[!is_gr1], ref_celltypes[ref_is_gr1])
  res_holdout1$comparison <- "holdout1_vs_obs1"
  res_holdout2 <- calculate_metrics(pred2, ref_obs2, celltypes[!is_gr2], ref_celltypes[ref_is_gr2])
  res_holdout2$comparison <- "holdout2_vs_obs2"
  list(res_train1, res_train2, res_holdout1, res_holdout2)  
}else{
  list(res_train1, res_train2)  
}

res %>%
  map(\(x) as.data.frame(x)) %>%
  bind_rows() %>%
  readr::write_tsv(file = outfile)



# Session Info
sessionInfo()




