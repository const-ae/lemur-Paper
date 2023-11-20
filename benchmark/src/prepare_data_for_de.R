library(SingleCellExperiment)
library(tidyverse)
source("src/utils/data_helper.R")
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Load data, subset it appropriately, and save it in a useful format")
pa <- argparser::add_argument(pa, "--dataset", type = "character", nargs = 1, help = "The name of a dataset") 
pa <- argparser::add_argument(pa, "--condition", type = "character", nargs = 1, help = "The name of condition") 
pa <- argparser::add_argument(pa, "--seed", type = "integer", nargs = 1, default = 1, help = "Seed to tame randomness") 
pa <- argparser::add_argument(pa, "--n_hvgs", type = "", nargs = "numeric", default = Inf, help = "The number of genes selected as highly variable") 
pa <- argparser::add_argument(pa, "--randomize", type = "character", nargs = 1, default = "cells", help = "Either randomize samples or cells") 
pa <- argparser::add_argument(pa, "--cut_at", type = "integer", nargs = Inf, help = "Determine how large the DE patches are") 
pa <- argparser::add_argument(pa, "--n_de_genes", type = "integer", nargs = 1, default = 100, help = "Determine how many genes are generated to have DE") 
pa <- argparser::add_argument(pa, "--lfc_mean", type = "numeric", nargs = Inf, help = "Determine how many genes are generated to have DE") 
pa <- argparser::add_argument(pa, "--clustering", type = "character", nargs = 1, default = "graph", help = "Determine how many genes are generated to have DE") 


pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --dataset angelidis
#                             --condition 3m
#                             --randomize cells
#                             --cut_at 1 10 25
#                             --n_de_genes 50
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

if(! pa$dataset %in% names(data_loader)){
  stop("No data loader defined for ", pa$dataset, ".\nThe options are ", toString(names(data_loader)))
}

randomization_levels <- match.arg(pa$randomization, c("cells", "samples"))
clustering_method <- match.arg(pa$clustering, c("graph", "kmeans", "celltype"))
if(all(is.na(pa$cut_at))){
  pa$cut_at <- c(1, 3, 10)
}
if(all(is.na(pa$lfc_mean))){
  pa$lfc_mean <- 1.5
}

output_file <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

config <- get_data_config(pa$dataset)
sce <- data_loader[[pa$dataset]]()

pa$condition <- paste0(pa$condition, collapse = " ")
sce <- sce[,colData(sce)[[config$main_covariate]] == pa$condition]

if(pa$n_hvgs < nrow(sce)){
  hvg <- order(-rowVars(assay(sce, config$assay_continuous)))
  n_hvgs <- min(nrow(sce), pa$n_hvgs)
  sce <- sce[hvg[seq_len(n_hvgs)],]
}

if(randomization_levels == "cells"){
  colData(sce)$fake_condition <- sample(c("fake_ctrl", "fake_trt"), size = ncol(sce), replace = TRUE)
  colData(sce)$sample <- colData(sce)[[config$sample_covariate]]
}else{
  sample_vec <- colData(sce)[[config$sample_covariate]]
  samples <- unique(sample_vec)
  trt_samples <- sample(samples, size = round(length(samples)), replace = FALSE)
  colData(sce)$fake_condition <- ifelse(sample_vec %in% trt_samples, "fake_trt", "fake_ctrl")
  colData(sce)$sample <- colData(sce)[[config$sample_covariate]]
}



# Generate DE
pca <- lemur:::pca(assay(sce, config$assay_continuous), n = min(nrow(sce), ncol(sce), 50))

if(clustering_method == "graph"){
  graph <- bluster::makeKNNGraph(t(pca$embedding), k = 15, BNPARAM = BiocNeighbors::AnnoyParam())
  clustering <- igraph::cluster_walktrap(graph)
}else if(clustering_method == "kmeans"){
  kmeans_clusterings <- lapply(unique(pa$cut_at), \(k){
    kmeans(t(pca$embedding), centers = k)$cluster
  })
  names(kmeans_clusterings) <- as.character(unique(pa$cut_at))
}else if(clustering_method == "celltype"){
  celltypes <- colData(sce)[[config$cell_type_column]]
  celltypes[is.na(celltypes)] <- "MISSING_ANNOTATION"
}


generate_gene_for_cluster <- function(n_clusters = 10, base_expr = -2, lfc_mean = 1, lfc_sd = 0.5, sample_sd = 0.1, overdispersion = 0.1, ...){
  cluster_assign <- if(clustering_method == "graph"){
    igraph::cut_at(clustering, no = n_clusters)
  }else if(clustering_method == "kmeans"){
    if(as.character(n_clusters) %in% names(kmeans_clusterings)){
      kmeans_clusterings[[as.character(n_clusters)]]
    }else{
      kmeans(t(pca$embedding), centers = n_clusters)$cluster
    }
  }else if(clustering_method == "celltype"){
    # Ignore 'n_clusters'
    celltypes
  }
  sel_cluster <- sample(unique(cluster_assign), size = 1)
  is_de_cell <- cluster_assign == sel_cluster
  
  eta_ctrl <- 0
  eta_trt <- rnorm(1, mean = lfc_mean, sd = lfc_sd)
  
  trt_eff <- colSums(lemur:::one_hot_encoding(sce$fake_condition)[c("fake_ctrl", "fake_trt"),] * c(eta_ctrl, eta_trt))
  mouse_mean <- rnorm(length(unique(sce$sample)), mean = 0, sd = sample_sd)
  mouse_eff <- colSums(lemur:::one_hot_encoding(sce$sample) * mouse_mean)
  sf <- colSums(counts(sce))
  sf <- sf / median(sf)
  mu <- 2^(trt_eff * is_de_cell + mouse_eff + base_expr) * sf
  counts <- rnbinom(n = ncol(sce), mu = mu, size = 1/overdispersion)
  list(n_clusters = n_clusters, sel_cluster = sel_cluster, is_de_cell = is_de_cell, 
       lfc = eta_trt - eta_ctrl, base_expr = base_expr,
       log_expression_level = unname(trt_eff * is_de_cell + mouse_eff + base_expr),
       mouse_mean = mouse_mean, counts = counts)
}


de_args <- tibble(name = paste0("simulated_gene-", seq_len(pa$n_de_genes)),
                  is_simulated =  TRUE,
                  cut_at = rep_len(pa$cut_at, length = pa$n_de_genes),
                  base_expr = runif(n = pa$n_de_genes, min = -7, max = 3), 
                  lfc = rep_len(pa$lfc_mean, length = pa$n_de_genes) * sample(c(-1, 1), size = pa$n_de_genes, replace = TRUE), 
                  is_de_cell = map(pa$n_de_genes, \(.) rep(FALSE, ncol(sce))))

new_counts <- matrix(0, nrow = pa$n_de_genes, ncol = ncol(sce))
for(idx in seq_len(pa$n_de_genes)){
  gv <- generate_gene_for_cluster(n_clusters = de_args$cut_at[idx], 
                                  base_expr = de_args$base_expr[idx],
                                  lfc_mean = de_args$lfc[idx],
                                  lfc_sd = 0, sample_sd = 0.1, overdispersion = 0.2)  
  new_counts[idx,] <- gv$counts
  de_args$is_de_cell[[idx]] <- gv$is_de_cell
}

neg_de_args <- tibble(name = rownames(sce), is_simulated = FALSE)

new_sce <- SingleCellExperiment(list(counts = rbind(unname(counts(sce)), new_counts)), colData = colData(sce),
                                rowData = bind_rows(neg_de_args, de_args))
colnames(new_sce) <- colnames(sce)
rownames(new_sce) <- rowData(new_sce)$name
assay(new_sce, "logcounts") <- transformGamPoi::shifted_log_transform(new_sce)


qs::qsave(new_sce, output_file)

# Session Info
sessionInfo()