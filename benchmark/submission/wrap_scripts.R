



prepare_data <- function(params, duration = "01:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/prepare_data.R", params = params, duration = duration, memory = memory)
}

#------------------------------------------------------------------------------------------------------

lemur_integration_prediction <- function(params, dep_job, duration = "03:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/lemur_integration_prediction.R", params = params,
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

linear_prediction <- function(params, dep_job, duration = "03:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/linear_prediction.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

harmony_integration <- function(params, dep_job, duration = "03:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/harmony_integration.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

pca_integration_prediction <- function(params, dep_job, duration = "03:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/pca_integration_prediction.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

cpa_integration_prediction <- function(params, dep_job, duration = "10:00:00", memory = "80GB"){
  MyWorkflowManager::wrap_script("src/run_cpa.py", params = params, 
                                 dependencies = list(dep_job), executor = "python",
                                 extra_args = "cpa_env", duration = duration, memory = memory)
}

#------------------------------------------------------------------------------------------------------

evaluate_integration <- function(params, dep_jobs, duration = "03:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/evaluate_integration.R", params = params, 
                                 dependencies = dep_jobs, duration = duration, memory = memory)
}

evaluate_prediction <- function(params, dep_jobs, duration = "10:00:00", memory = "80GB"){
  MyWorkflowManager::wrap_script("src/evaluate_prediction.R", params = params, 
                                 dependencies = dep_jobs, duration = duration, memory = memory)
}






