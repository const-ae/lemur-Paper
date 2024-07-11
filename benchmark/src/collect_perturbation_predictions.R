library(SingleCellExperiment)
library(tidyverse)

pa <- argparser::arg_parser("Collect perturbation predictions")
pa <- argparser::add_argument(pa, "--job_ids", type = "character", nargs = Inf, help = "The job ids") 
pa <- argparser::add_argument(pa, "--names", type = "character", nargs = Inf, help = "The method names") 

pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --names norman-1-scgpt norman-1-scgpt_oneshot  norman-1-gears norman-1-ground_truth norman-1-additive_model norman-1-pylemur
#                             --job_ids 7f3516ac013be-2e77bf6429630 7f3516ac013be-bffe4cf3eea4c b3584d8d39033-2e77bf6429630 9c0ba4ceae23e-2e77bf6429630 6d8f405477fb2-2e77bf6429630 bab2cbb8f6965-2e77bf6429630
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)

out_dir <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

stopifnot(length(pa$names) == length(pa$job_ids))

safe_as_numeric = function(x, na.strings = "NA") {
  if(is.character(x)){
    na = x %in% na.strings
    x[na] = "0"
    x = as.numeric(x)
    x[na] = NA_real_
    x
  }else{
    as.numeric(x)
  }
}

raw_res <- map2(pa$job_ids, pa$names, \(id, name){
    res <- rjson::fromJSON(file = file.path(pa$working_dir, "results", id, "all_predictions.json"))
    params_id <- stringr::str_split_fixed(id, "-", n = 2)[,2]
    params <- yaml::read_yaml(file.path(pa$working_dir, "params", params_id))
    print(paste0("Finished reading ", name))
    
    list(id = id, name = name, perturbation = names(res), prediction = map(res, safe_as_numeric))
  }) 

parameters <- map2(pa$job_ids, pa$names, \(id, name){
  params_id <- stringr::str_split_fixed(id, "-", n = 2)[,2]
  params <- yaml::read_yaml(file.path(pa$working_dir, "params", params_id))
  test_train_labels <- rjson::fromJSON(file = file.path(pa$working_dir, "results", params$test_train_config_id))
  print(paste0("Finished reading ", name))
  
  list(id = id, name = name, parameters = params, test_train_labels = test_train_labels)
})


# Store output
tdir <- file.path(tempdir(), paste0("temp-", pa$result_id))
dir.create(tdir)
readr::write_rds(raw_res, file.path(tdir, "predictions.RDS"), compress = "gz")
readr::write_rds(parameters, file.path(tdir, "parameters.RDS"), compress = "gz")
file.rename(tdir, out_dir)

#### Session Info
sessionInfo()


