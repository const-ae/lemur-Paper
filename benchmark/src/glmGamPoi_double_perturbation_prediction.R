library(SingleCellExperiment)
library(tidyverse)
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Run glmGamPoi for double perturbation prediction")
pa <- argparser::add_argument(pa, "--dataset_name", type = "character", help = "") 
pa <- argparser::add_argument(pa, "--test_train_config_id", type = "character", nargs = 1, help = "") 
pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The seed")
pa <- argparser::add_argument(pa, "--ridge_penalty", type = "numeric", default = 0, nargs = 1, help = "The ridge penalty on non-ctrl columns")


pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --dataset_name norman
#                             --test_train_config_id 8443ed21d2ac4-f8716281f960b
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

out_dir <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

sce <- zellkonverter::readH5AD(file.path("data/gears_pert_data/", pa$dataset_name, "perturb_processed.h5ad"), reader = "R")
colData(sce) <- colData(sce) %>%
  as_tibble() %>%
  mutate(pert_split = str_split(condition, "[+_]")) %>%
  mutate(pert_split = map(pert_split, \(x){
    if(all(x == "ctrl") || all(x == "")) "ctrl"
    else if(length(x) == 1) sort(c(x, "ctrl"))
    else sort(x)
  })) %>%
  mutate(pert = map_chr(pert_split, paste0, collapse = "+")) %>%
  as.data.frame() %>%
  DataFrame()


train_info <- rjson::fromJSON(file = file.path(pa$working_dir, "results", pa$test_train_config_id))
train_df <- tibble(training = rep(names(train_info), lengths(train_info)), perturbation = unlist(train_info)) %>%
  mutate(pert_split = str_split(perturbation, "[+_]")) %>%
  mutate(pert_split = map(pert_split, \(x){
    if(all(x == "ctrl") || all(x == "")) "ctrl"
    else if(length(x) == 1) sort(c(x, "ctrl"))
    else sort(x)
  })) %>%
  mutate(pert = map_chr(pert_split, paste0, collapse = "+"))

train_sce <- sce[,sce$pert %in% filter(train_df, training == "train")$pert]


se <- glmGamPoi::pseudobulk(train_sce, group_by = vars(pert))
make_design <- function(conds, levels = NULL){
  if(is.null(levels)){
    levels <- sort(unique(unlist(conds)))
  }
  mat <- matrix(0, nrow = length(conds), ncol = length(levels))
  colnames(mat) <- levels
  for(idx in seq_along(conds)){
    mat[idx, levels %in% conds[[idx]]] <- 1
  }
  mat
}

des <- make_design(str_split(se$pert,  "\\+"))
des[,"ctrl"] <- 1

ridge_penalty <- pmax(rep(pa$ridge_penalty, ncol(des)), 1e-6)
ridge_penalty[which(colnames(des) == "ctrl")] <- 0

attr(des, "ignore_degeneracy") <- TRUE
fit <- glmGamPoi::glm_gp(se, design = des, size_factors = "ratio", verbose = TRUE, 
                         overdispersion_shrinkage = FALSE, overdispersion = 0.01, 
                         ridge_penalty = ridge_penalty)

uniq_pert <- unique(colData(sce)$pert)
preds <- lapply(str_split(uniq_pert, "\\+"), \(conds){
  drop(predict(fit, newdata = make_design(list(c(conds, "ctrl")), levels = colnames(des)), type = "link"))
})
names(preds) <- uniq_pert

tmp_out_dir <- paste0(out_dir, "-tmp")
dir.create(tmp_out_dir)
readr::write_file(rjson::toJSON(preds, indent = 4), file.path(tmp_out_dir, "all_predictions.json"))
file.rename(tmp_out_dir, out_dir)

#### Session Info
sessionInfo()


