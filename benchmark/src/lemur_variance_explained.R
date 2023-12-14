library(tidyverse)
library(SingleCellExperiment)
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Analyze fraction of variance explained as a function of latent dimensions")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--dataset_config", type = "character", nargs = 1, help = "The name of a dataset in datasets.yaml") 
pa <- argparser::add_argument(pa, "--dataset_config_override", type = "character", default = "", nargs = Inf, help = "Override settings from datasets.yaml") 

pa <- argparser::add_argument(pa, "--pca_dims", type = "integer", nargs = Inf, help = "The number of PCA dimensions")
pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The number of PCA dimensions")

pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --data_id f1809110732c2-e5f290978d2e2
#                             --dataset_config kang
#                             --pca_dims 3 10 30
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

config <- get_data_config(pa$dataset_config, pa$dataset_config_override)

out_file <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

sce <- zellkonverter::readH5AD(file.path(pa$working_dir, "results", pa$data_id, "train.h5ad"))

mat <- assay(sce, config$assay_continuous)
message("Total RSS: ", sprintf("%.3g", sum(mat^2)))
resid <- ResidualMatrix::ResidualMatrix(t(mat), design = model.matrix(as.formula(paste0("~ ", config$main_covariate)), data = colData(sce)))
mat <- t(as.matrix(resid))
message("Total RSS (centered): ", sprintf("%.3g", sum(mat^2)))

var_explained <- function(Y, obj, FUN, subset = NULL, ...){
  if(is.null(subset)){
    1 - sum((Y - FUN(obj, ...))^2) / sum(Y^2)
  }else{
    1 - sum((Y[,subset,drop=FALSE] - FUN(obj, subset = subset, ...))^2) / sum(Y[,subset,drop=FALSE]^2)
  }
}

predict_lemur <- function(fit, subset = NULL){
  if(is.null(subset)){
    predict(fit)
  }else{
    predict(fit[,subset])
  }
}
predict_pca <- function(fit, subset = NULL){
  if(is.null(subset)){
    fit$rotation %*% t(fit$x)
  }else{
    fit$rotation %*% t(fit$x[subset,,drop=FALSE])
  }
}

pred_func <- list(pca = predict_pca, lemur = predict_lemur)

subsets <- c(list(NULL), lapply(config$contrast, \(lvl) colData(sce)[[config$main_covariate]] == lvl))
names(subsets) <- c("all", config$contrast)

dataset_overview <- tibble(dataset = pa$dataset_config, n_genes = nrow(mat), n_cells = ncol(mat))

res <- bind_rows(lapply(pa$pca_dims, \(dim){
  message("Starting: ", dim)
  
  pca_time <- system.time({
    pca_fit <- irlba::prcomp_irlba(t(mat), n = dim, center = FALSE, scale = FALSE)
  })
  
  lemur_time <- system.time({
    lemur_fit <- lemur::lemur(mat, design = as.formula(paste0("~ ", config$main_covariate)), col_data = colData(sce),
                 n_embedding = dim, test_fraction = 0, linear_coefficient_estimator = "zero", verbose = FALSE)
  })
  lemur_al_landmark_time <- system.time({
    lemur_fit_al1 <- lemur::align_by_grouping(lemur_fit, grouping = colData(sce)[[config$cell_type_column]])
  })
  lemur_al_harmony_time <- tryCatch({
    system.time({
      lemur_fit_al2 <- lemur::align_harmony(lemur_fit)
    })
    },
    error = function(err){
      lemur_al_landmark_time * NA
    })
  
  models <- list(pca = pca_fit, lemur = lemur_fit)
  
  timing_df <- bind_rows(as.list(pca_time), as.list(lemur_time), 
                         as.list(lemur_time + lemur_al_landmark_time), 
                         as.list(lemur_time + lemur_al_harmony_time)) %>%
    mutate(method = c("pca", "lemur", "lemur_landmark", "lemur_harmony"),
           dimensions = dim)
  
  tidyr::expand_grid(dimensions = dim, 
                     method = c("pca", "lemur"),
                     subset = c("all", config$contrast)) %>%
    mutate(var_expl = map2_dbl(method, subset, \(m, s){
      var_explained(mat, models[[m]], pred_func[[m]], subset = subsets[[s]])
    }))  %>%
    full_join(timing_df, by = c("method", "dimensions")) %>%
    bind_cols(dataset_overview)
}))

write_tsv(res, out_file)


#### Session Info
sessionInfo()


