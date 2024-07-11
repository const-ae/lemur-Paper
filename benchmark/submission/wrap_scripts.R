



prepare_data <- function(params, duration = "06:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/prepare_data.R", params = params, duration = duration, memory = memory)
}

prepare_de_data <- function(params, duration = "06:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/prepare_data_for_de.R", params = params, duration = duration, memory = memory)
}

prepare_perturbation_data <- function(params, duration = "01:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/prepare_perturbation_data.py", params = params, executor = "python", 
                                 extra_args = "gears_env2", duration = duration, memory = memory)
}


#------------------------------------------------------------------------------------------------------

lemur_integration_prediction <- function(params, dep_job, duration = "06:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/lemur_integration_prediction.R", params = params,
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

linear_prediction <- function(params, dep_job, duration = "06:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/linear_prediction.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

no_change_prediction <- function(params, dep_job, duration = "02:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/no_change_prediction.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

harmony_integration <- function(params, dep_job, duration = "06:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/harmony_integration.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

pca_integration_prediction <- function(params, dep_job, duration = "06:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/pca_integration_prediction.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

cpa_on_counts_integration_prediction <- function(params, dep_job, duration = "24:00:00", memory = "80GB"){
  MyWorkflowManager::wrap_script("src/run_cpa_on_counts.py", params = params, 
                                 dependencies = list(dep_job), executor = "python",
                                 extra_args = "cpa_env", duration = duration, memory = memory)
}

scvi_integration_prediction <- function(params, dep_job, duration = "24:00:00", memory = "80GB"){
  MyWorkflowManager::wrap_script("src/run_scvi.py", params = params, 
                                 dependencies = list(dep_job), executor = "python",
                                 extra_args = "scvi_env2", duration = duration, memory = memory)
}

#------------------------------------------------------------------------------------------------------


lemur_de <- function(params, dep_job, duration = "06:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/lemur_differential_expression.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

global_de <- function(params, dep_job, duration = "01:00:00", memory = "20GB"){
  MyWorkflowManager::wrap_script("src/global_differential_expression.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

celltype_de <- function(params, dep_job, duration = "01:00:00", memory = "20GB"){
  MyWorkflowManager::wrap_script("src/celltype_differential_expression.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

cluster_de <- function(params, dep_job, duration = "01:00:00", memory = "20GB"){
  MyWorkflowManager::wrap_script("src/cluster_differential_expression.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}

milo_de <- function(params, dep_job, duration = "01:00:00", memory = "20GB"){
  MyWorkflowManager::wrap_script("src/miloDE_differential_expression.R", params = params, 
                                 dependencies = list(dep_job),
                                 duration = duration, memory = memory)
}


evaluate_de_results <- function(params, dep_jobs, duration = "05:00:00", memory = "60GB"){
  MyWorkflowManager::wrap_script("src/evaluate_de_results.R", params = params, 
                                 dependencies = dep_jobs,
                                 duration = duration, memory = memory)
}


#------------------------------------------------------------------------------------------------------

evaluate_integration <- function(params, dep_jobs, duration = "06:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/evaluate_integration.R", params = params, 
                                 dependencies = dep_jobs, duration = duration, memory = memory)
}

evaluate_prediction <- function(params, dep_jobs, duration = "24:00:00", memory = "80GB"){
  MyWorkflowManager::wrap_script("src/evaluate_prediction.R", params = params, 
                                 dependencies = dep_jobs, duration = duration, memory = memory)
}


visualize_integration <- function(params, dep_jobs, duration = "03:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/visualize_integration.R", params = params, 
                                 dependencies = dep_jobs, duration = duration, memory = memory)
}


#------------------------------------------------------------------------------------------------------

calc_variance_explained <- function(params, dep_jobs, duration = "05:00:00", memory = "60GB"){
  MyWorkflowManager::wrap_script("src/lemur_variance_explained.R", params = params, 
                                 dependencies = dep_jobs,
                                 duration = duration, memory = memory)
}


#------------------------------------------------------------------------------------------------------

scgpt_combinatorial_prediction <- function(params, dep_jobs, duration = "10:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/run_scgpt.py", params = params, 
                                 dependencies = dep_jobs, executor = "python", 
                                 extra_args = "flashattn_env",
                                 extra_slurm_arg = "-C gaming  -p gpu-el8 --gpus 1 --ntasks-per-gpu 1",
                                 duration = duration, memory = memory)
}

gears_combinatorial_prediction <- function(params, dep_jobs, duration = "10:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/run_gears.py", params = params, 
                                 dependencies = dep_jobs, executor = "python", 
                                 extra_args = "gears_env2",
                                 extra_slurm_arg = "-C gaming  -p gpu-el8 --gpus 1 --ntasks-per-gpu 1",
                                 duration = duration, memory = memory)
}

additive_model_combinatorial_prediction <- function(params, dep_jobs, duration = "01:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/run_additive_model.py", params = params, 
                                 dependencies = dep_jobs, executor = "python", 
                                 extra_args = "gears_env2",
                                 duration = duration, memory = memory)
}

glmGamPoi_combinatorial_prediction <- function(params, dep_jobs, duration = "03:00:00", memory = "60GB"){
  MyWorkflowManager::wrap_script("src/glmGamPoi_double_perturbation_prediction.R", params = params, 
                                 dependencies = dep_jobs, executor = "R", 
                                 duration = duration, memory = memory)
}

pylemur_combinatorial_prediction <- function(params, dep_jobs, duration = "03:00:00", memory = "60GB"){
  MyWorkflowManager::wrap_script("src/pylemur_double_perturbation_prediction.py", params = params, 
                                 dependencies = dep_jobs, executor = "python", 
                                 extra_args = "pylemur_env",
                                 duration = duration, memory = memory)
}

ground_truth_combinatorial_prediction <- function(params, dep_jobs, duration = "01:00:00", memory = "40GB"){
  MyWorkflowManager::wrap_script("src/run_ground_truth_for_combinatorial_perturbations.py", params = params, 
                                 dependencies = dep_jobs, executor = "python", 
                                 extra_args = "gears_env2",
                                 duration = duration, memory = memory)
}


linear_perturbation_prediction <- function(params, dep_jobs, duration = "03:00:00", memory = "60GB"){
  MyWorkflowManager::wrap_script("src/run_linear_perturbation_prediction.R", params = params, 
                                 dependencies = dep_jobs, executor = "R", 
                                 duration = duration, memory = memory)
}

transfer_perturbation_prediction <- function(params, dep_jobs, duration = "03:00:00", memory = "60GB"){
  MyWorkflowManager::wrap_script("src/run_transfer_perturbation_prediction.R", params = params, 
                                 dependencies = dep_jobs, executor = "R", 
                                 duration = duration, memory = memory)
}

collect_perturbation_predictions <- function(params, dep_jobs, duration = "03:00:00", memory = "60GB"){
  MyWorkflowManager::wrap_script("src/collect_perturbation_predictions.R", params = params, 
                                 dependencies = dep_jobs, executor = "R", 
                                 duration = duration, memory = memory)
}


