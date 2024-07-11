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


mmd_distance_wrapper <- function(data, ref_data, is_gr1, ref_is_gr1, gammas = 10^seq(1, -3, length.out = 50), n_cells = c(100, 200, 500)){
  stopifnot(length(is_gr1) == ncol(data))
  stopifnot(length(ref_is_gr1) == ncol(ref_data))
  stopifnot(is.logical(is_gr1) && is.logical(ref_is_gr1))
  ex <- data[, is_gr1,drop=FALSE]
  ey <- data[,!is_gr1,drop=FALSE]
  rx <- ref_data[, ref_is_gr1,drop=FALSE]
  ry <- ref_data[,!ref_is_gr1,drop=FALSE]
  
  res <- rep(NA, length(gammas) * length(n_cells))
  idx <- 1
  
  for(nc in n_cells){
    if(nc > min(ncol(ex), ncol(ey), ncol(rx), ncol(ry))){
      break
    }
    sel_ex <- sample.int(ncol(ex), size = nc)
    sel_ey <- sample.int(ncol(ey), size = nc)
    sel_rx <- sample.int(ncol(rx), size = nc)
    sel_ry <- sample.int(ncol(ry), size = nc)
    ext_sub <- t(ex[,sel_ex,drop=FALSE])
    eyt_sub <- t(ey[,sel_ey,drop=FALSE])
    rxt_sub <- t(rx[,sel_rx,drop=FALSE])
    ryt_sub <- t(ry[,sel_ry,drop=FALSE])
    for(g in gammas){
      res[idx] <- (mmd_distance(ext_sub, ryt_sub, gamma = g) + mmd_distance(rxt_sub, eyt_sub, gamma = g)) / 2
      idx <- idx + 1
    }
  }
  mean(res, na.rm = TRUE)
}


wasserstein_wrapper <- function(data, ref_data, is_gr1, ref_is_gr1, n_cells = c(100, 200, 500)){
  stopifnot(length(is_gr1) == ncol(data))
  stopifnot(length(ref_is_gr1) == ncol(ref_data))
  stopifnot(is.logical(is_gr1) && is.logical(ref_is_gr1))
  
  ex <- data[, is_gr1,drop=FALSE]
  ey <- data[,!is_gr1,drop=FALSE]
  rx <- ref_data[, ref_is_gr1,drop=FALSE]
  ry <- ref_data[,!ref_is_gr1,drop=FALSE]
  
  res <- rep(NA, length(n_cells))
  idx <- 1
  
  for(nc in n_cells){
    if(nc > min(ncol(ex), ncol(ey), ncol(rx), ncol(ry))){
      break
    }
    sel_ex <- sample.int(ncol(ex), size = nc)
    sel_ey <- sample.int(ncol(ey), size = nc)
    sel_rx <- sample.int(ncol(rx), size = nc)
    sel_ry <- sample.int(ncol(ry), size = nc)
    ext_sub <- t(ex[,sel_ex,drop=FALSE])
    eyt_sub <- t(ey[,sel_ey,drop=FALSE])
    rxt_sub <- t(rx[,sel_rx,drop=FALSE])
    ryt_sub <- t(ry[,sel_ry,drop=FALSE])
    
    tryCatch({
      res[idx] <- (transport::wasserstein(transport::pp(ext_sub), transport::pp(ryt_sub)) +
                     transport::wasserstein(transport::pp(rxt_sub), transport::pp(eyt_sub))) / 2
    }, error = function(err){
      message("Optimal transport failed for ", nc)
    })
    idx <- idx + 1
  }
  mean(res, na.rm = TRUE)
}

mixing_score <- function(data, ref_data, is_gr1, ref_is_gr1, k = 20){
  stopifnot(length(is_gr1) == ncol(data))
  stopifnot(length(ref_is_gr1) == ncol(ref_data))
  stopifnot(is.logical(is_gr1) && is.logical(ref_is_gr1))
  stopifnot(length(k) == 1 && k > 0 && k < ncol(ref_data))
  
  ref_index <- BiocNeighbors::buildAnnoy(ref_data, transposed = TRUE)
  neigh <- BiocNeighbors::queryAnnoy(precomputed = ref_index, query = data, k = k, transposed = TRUE, 
                                     get.index = TRUE, get.distance = FALSE)$index
  
  mean(vapply(seq_len(ncol(data)), \(idx){
    true_label <- is_gr1[idx]
    sum(true_label != ref_is_gr1[neigh[idx,]])
  }, FUN.VALUE = numeric(1L)))
}

overlap_score <- function(data, noint_data, is_gr1, k = 20){
  score1 <- overlap_score_impl(ref_embedding[,is_gr1], noint_embedding[,is_gr1])
  score2 <- overlap_score_impl(ref_embedding[,!is_gr1], noint_embedding[,!is_gr1])
  list(overlap = (score1 + score2) / 2)
}

overlap_score_impl <- function(data, noint_data, k = 20){
  stopifnot(ncol(noint_data) == ncol(data))
  stopifnot(length(k) == 1 && k > 0 && k < ncol(data))
  
  noint_index <- BiocNeighbors::buildAnnoy(noint_data, transposed = TRUE)
  noint_nei <- BiocNeighbors::queryAnnoy(precomputed = noint_index, query = noint_data, k = k, transposed = TRUE, 
                                         get.index = TRUE, get.distance = FALSE)$index
  data_index <- BiocNeighbors::buildAnnoy(data, transposed = TRUE)
  data_nei <- BiocNeighbors::queryAnnoy(precomputed = data_index, query = data, k = k, transposed = TRUE, 
                                          get.index = TRUE, get.distance = FALSE)$index
  
  mean(vapply(seq_len(ncol(data)), \(idx){
    length(intersect(noint_nei[idx,], data_nei[idx,]))
  }, FUN.VALUE = numeric(1L)))
}


structure_conservation <- function(large_data, embedding, is_gr1){
  pca_ref_gr1 <- irlba::prcomp_irlba(t(large_data[, is_gr1]), n = 30)  
  clustering_ref_gr1 <- bluster::clusterRows(pca_ref_gr1$x, bluster::KNNGraphParam())

  clustering_emb <- bluster::clusterRows(t(embedding[,is_gr1]), bluster::KNNGraphParam(), 
                                         full = TRUE)
  clustering_emb_cut <- igraph::cut_at(clustering_emb$objects$communities, 
                                       no = length(unique(clustering_ref_gr1)))
  list(
    ARI_onlog = aricode::ARI(clustering_ref_gr1, clustering_emb_cut),
    AMI_onlog = aricode::AMI(clustering_ref_gr1, clustering_emb_cut),
    NMI_onlog = aricode::NMI(clustering_ref_gr1, clustering_emb_cut)
  )
}

structure_conservation2 <- function(noint_embedding, embedding, is_gr1){
  clustering_ref_gr1 <- bluster::clusterRows(t(noint_embedding[,is_gr1]), bluster::KNNGraphParam())
  clustering_emb <- bluster::clusterRows(t(embedding[,is_gr1]), bluster::KNNGraphParam(), 
                                         full = TRUE)
  clustering_emb_cut <- igraph::cut_at(clustering_emb$objects$communities, 
                                       no = length(unique(clustering_ref_gr1)))
  list(
    ARI_onself = aricode::ARI(clustering_ref_gr1, clustering_emb_cut),
    AMI_onself = aricode::AMI(clustering_ref_gr1, clustering_emb_cut),
    NMI_onself = aricode::NMI(clustering_ref_gr1, clustering_emb_cut)
  )
}


signal_to_noise_ratio <- function(embedding, is_gr1, cell_types){
  overall_ss <- sum(residuals(lm(t(embedding) ~ 1))^2)
  celltype_ss <- sum(residuals(lm(t(embedding) ~ cell_types))^2)
  gr1_plus_celltype_ss <- sum(residuals(lm(t(embedding) ~ cell_types + is_gr1))^2)
  celltype_ss / overall_ss
  (celltype_ss - gr1_plus_celltype_ss) / overall_ss
  
  list(celltypeSS_over_totalSS = celltype_ss / overall_ss, batchSS_over_totalSS = (celltype_ss - gr1_plus_celltype_ss) / overall_ss)
}

calculate_metrics <- function(data, ref_data, is_gr1, ref_is_gr1){
  purrr::flatten(list(
    mmd = mmd_distance_wrapper(data, ref_data, is_gr1, ref_is_gr1),
    wasserstein = wasserstein_wrapper(data, ref_data, is_gr1, ref_is_gr1),
    mix = mixing_score(data, ref_data, is_gr1, ref_is_gr1)
  ))
}


pa <- argparser::arg_parser("Evaluate integration results")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--dataset_config", type = "character", nargs = 1, help = "The name of a dataset in datasets.yaml") 
pa <- argparser::add_argument(pa, "--dataset_config_override", type = "character", default = "", nargs = Inf, help = "Override settings from datasets.yaml") 
pa <- argparser::add_argument(pa, "--integration_id", type = "character", help = "The result id of the integration job")
pa <- argparser::add_argument(pa, "--reference_id", type = "character", help = "The result id of the reference job")

pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The seed")

pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --data_id 337dc3137fe9d-e5f290978d2e2
#                             --dataset_config kang
#                             --integration_id 7a118bf7498c0-33b38b15c9f6e
#                             --reference_id 42a82ea22a4a9-33b38b15c9f6e
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

ref_is_gr1 <- colData(sce)[[config$main_covariate]] == config$contrast[1]
ref_celltypes <- colData(sce)[[config$cell_type_column]]
ref_embedding <- as.matrix(data.table::fread(train_file, sep = "\t", header = FALSE, col.names = colnames(sce)))

noint_embedding <- file.path(pa$working_dir, "results", pa$reference_id, "train-embedding.tsv") |>
  data.table::fread(sep = "\t", header = FALSE, col.names = colnames(sce)) |>
  as.matrix()

# Bring all data to the same scale
rm <- rowMeans(ref_embedding)
ref_embedding <- ref_embedding - rm
scale_factor <- mean(sqrt(colSums(ref_embedding^2)))
ref_embedding <- ref_embedding / scale_factor + rm

res1 <- calculate_metrics(ref_embedding, ref_embedding, ref_is_gr1, ref_is_gr1)
res1 <- c(res1, structure_conservation(assay(sce, config$assay_continuous), ref_embedding, ref_is_gr1))
res1 <- c(res1, structure_conservation2(noint_embedding, ref_embedding, ref_is_gr1))
res1 <- c(res1, overlap_score(ref_embedding, noint_embedding, ref_is_gr1))
res1 <- c(res1, signal_to_noise_ratio(ref_embedding, ref_is_gr1, ref_celltypes))
res1$comparison <- "train_vs_train"

res <- if(file.exists(holdout_file)){
  embedding <- as.matrix(data.table::fread(holdout_file, sep = "\t", header = FALSE, col.names = colnames(holdout)))
  embedding <- (embedding - rm) / scale_factor + rm
  noint_holdout_embedding <- file.path(pa$working_dir, "results", pa$reference_id, "holdout-embedding.tsv") |>
    data.table::fread(sep = "\t", header = FALSE, col.names = colnames(holdout)) |>
    as.matrix()
  
  is_gr1 <- colData(holdout)[[config$main_covariate]] == config$contrast[1]
  celltypes <- colData(holdout)[[config$cell_type_column]]
  res2 <- calculate_metrics(embedding, ref_embedding, is_gr1, ref_is_gr1)  
  res2 <- c(res2, structure_conservation(assay(holdout, config$assay_continuous), embedding, is_gr1))
  res2 <- c(res2, structure_conservation2(noint_holdout_embedding, embedding, is_gr1))
  res2 <- c(res2, overlap_score(embedding, noint_holdout_embedding, is_gr1))
  res2 <- c(res2, signal_to_noise_ratio(embedding, is_gr1, celltypes))
  
  res2$comparison <- "holdout_vs_train"  
  rbind(as.data.frame(res1), 
        as.data.frame(res2))
}else{
  as.data.frame(res1)
}


readr::write_tsv(res, file = outfile)



# Session Info
sessionInfo()
