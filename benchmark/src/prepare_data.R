source("src/utils/data_helper.R")

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

valid_variations <- c("random_holdout", "identity")
var <- pa$variation
if(! var %in% valid_variations){
  stop("Variation ", var,  " is not a recognized dataset variation")
}

output_file <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------



sce <- data_loader[[pa$dataset]]()
if(var == "identity"){
  train <- sce
  holdout <- NULL
}else if(var == "random_holdout"){
  holdout_lgl <- runif(seq_len(ncol(sce))) < pa$hold_out_percentage
  train <- sce[,!holdout_lgl]
  holdout <- sce[,holdout_lgl]
}

# Save everything
qs::qsave(list(train = train, holdout = holdout), output_file)

# Session Info
sessionInfo()