source("src/utils/data_helper.R")
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Run pseudobulk differential expression per cluster")
pa <- argparser::add_argument(pa, "--dataset", type = "character", nargs = 1, help = "The name of a dataset") 
pa <- argparser::add_argument(pa, "--variation", type = "character", nargs = 1, help = "The data variation that is generated") 
pa <- argparser::add_argument(pa, "--seed", type = "integer", nargs = 1, default = 1, help = "Seed to tame randomness") 
pa <- argparser::add_argument(pa, "--hold_out_percentage", type = "", nargs = "numeric", default = 0.2, help = "The holdout percentage if variation is 'random_holdout'") 
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

valid_variations <- c("random_holdout", "identity", "small")
var <- pa$variation
if(! var %in% valid_variations){
  stop("Variation ", var,  " is not a recognized dataset variation")
}

output_file <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

config <- get_data_config(pa$dataset)


sce <- data_loader[[pa$dataset]]()
if(var == "identity"){
  train <- sce
  holdout <- NULL
}else if(var == "random_holdout"){
  holdout_lgl <- runif(seq_len(ncol(sce))) < pa$hold_out_percentage
  train <- sce[,!holdout_lgl]
  holdout <- sce[,holdout_lgl]
}else if(var == "small"){
  hvg <- order(-rowVars(assay(sce, config$assay_continuous)))
  n_genes <- min(nrow(sce), 500)
  n_cells <- min(ncol(sce), 5000)
  train <- sce[hvg[seq_len(n_genes)], sample.int(ncol(sce), size = n_cells)]
  holdout <- sce[hvg[seq_len(n_genes)], sample.int(ncol(sce), size = n_cells)]
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