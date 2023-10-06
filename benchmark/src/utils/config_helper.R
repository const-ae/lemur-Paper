

get_data_config <- function(dataset_name = NULL, overrides = NULL){
  config <- yaml::read_yaml("data/datasets.yaml")
  if(is.null(dataset_name) || all(dataset_name == "")){
    # do nothing
  }else{
    if(! dataset_name %in% names(config)){
      stop("No data config defined for ", dataset_name, ".\nThe options are ", toString(names(config)))
    }
    config <- config[[dataset_name]]
  }
  
  if(is.null(overrides) || all(overrides == "")){
    config
  }else{
    update_config(config, overrides)
  }
}



update_config <- function(config, update){
  lhs_rhs <- stringr::str_split_fixed(update, "=", n = 2)
  for(idx in seq_len(nrow(lhs_rhs))){
    lhs_elem <- as.list(stringr::str_split_1(lhs_rhs[idx,1], "\\."))
    purrr::pluck(config, !!!lhs_elem) <- yaml::read_yaml(text = lhs_rhs[idx,2])
  }
  config
}
