

lemur_script <- function(params, duration = "01:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("scripts/de_lemur.R", params = params, duration = duration, memory = memory)
}


prepare_data <- function(params, duration = "01:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/prepare_data.R", params = params, duration = duration, memory = memory)
}


linear_prediction <- function(params, dep_job, duration = "01:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/linear_prediction.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

harmony_integration <- function(params, dep_job, duration = "01:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/harmony_integration.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

pca_integration_prediction <- function(params, dep_job, duration = "01:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/pca_integration_prediction.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

cpa_integration_prediction <- function(params, dep_job, duration = "01:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/run_cpa.py", params = params, 
                                 dependencies = list(dep_job), executor = "python",
                                 extra_args = "cpa_env", duration = duration, memory = memory)
}
