import os
import sys
print(os.getcwd())
sys.path.insert(0, os.getcwd())  


import argparse
import src.utils.config_helper as ch
from pathlib import Path


# Import AnnData object with scanpy
import scanpy as sc
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
import cpa
import session_info

print("Hello from python")

parser = argparse.ArgumentParser(description='Run CPA')
parser.add_argument('--data_id', dest='data_id', action='store', required = True, help='The id of a file in output/results')
parser.add_argument("--dataset_config", dest = "dataset_config", required = True, help = "The name of a dataset in datasets.yaml")
parser.add_argument("--dataset_config_override", dest = "dataset_config_override", default = [], nargs = "*", help = "Override settings from datasets.yaml")

parser.add_argument('--n_latent', dest='n_latent', action='store', default = 15, type = int,
                    help='The number of latent dimensions')
parser.add_argument('--max_epochs', dest='max_epochs', action='store', default = 100, type = int,
                    help='The number of training epochs.')

parser.add_argument("--working_dir", dest = "working_dir", action='store', required = True, help = "The directory that contains the params, results, scripts etc.")
parser.add_argument("--result_id", dest = "result_id", action='store', required = True, help = "The result_id")
args = parser.parse_args()
# args = parser.parse_args(["--data_id", "8f872edd6a7eb-827531ae33574",
#     "--dataset_config", "kang", "--n_latent", "20", "--working_dir",
#     "/scratch/ahlmanne/lemur_benchmark", "--result_id", "0"])
print(args)


config = ch.get_data_config(args.dataset_config, args.dataset_config_override)
out_dir = args.working_dir + "/results/" + args.result_id

# --------------------------------------------------------

adata = sc.read_h5ad(args.working_dir + "/results/" +  args.data_id + "/train.h5ad")
holdout_adata = sc.read_h5ad(args.working_dir + "/results/" +  args.data_id + "/holdout.h5ad")

if is_numeric_dtype(adata.obs[config['main_covariate']]):
  adata.obs['key'] = np.where(adata.obs[config['main_covariate']] == config['contrast'][0], "ctrl", "pert")
  holdout_adata.obs['key'] = np.where(holdout_adata.obs[config['main_covariate']] == config['contrast'][0], "ctrl", "pert")
  config['num_covariate'] = config['main_covariate']
  config['main_covariate'] = 'key'
  config['num_contrast'] = config['contrast']
  config['contrast'] = ["ctrl", "pert"]

  cpa.CPA.setup_anndata(adata, perturbation_key=config['main_covariate'],
       control_group = config['contrast'][0], is_count_data = False, dosage_key = config['num_covariate'])
  cpa.CPA.setup_anndata(holdout_adata, perturbation_key=config['main_covariate'],
       control_group = config['contrast'][0], is_count_data = False, dosage_key = config['num_covariate'])
else:
  cpa.CPA.setup_anndata(adata, perturbation_key=config['main_covariate'],
       control_group = config['contrast'][0], is_count_data = False)
  cpa.CPA.setup_anndata(holdout_adata, perturbation_key=config['main_covariate'],
       control_group = config['contrast'][0], is_count_data = False)


model = cpa.CPA(adata, n_latent = args.n_latent, recon_loss = 'gauss')
model.train(max_epochs= args.max_epochs)

# Get integrated embeddings
latent_outputs = model.get_latent_representation(adata, batch_size=1024)
basal_predict = model.custom_predict(basal = True, add_batch = False, add_pert = False, adata = holdout_adata)


# Make predictions
conditions = adata.obs[config['main_covariate']].unique()
trained_preds = []
holdout_preds = []
for cond in conditions:
  if cond == config['contrast'][0]:
    train_pred = model.custom_predict(basal = True, add_batch = False, add_pert = False)
    holdout_pred = model.custom_predict(basal = True, add_batch = False, add_pert = False, adata = holdout_adata)
  else:
    train_pred = model.custom_predict(basal = False, add_batch = False, add_pert = True)
    holdout_pred = model.custom_predict(basal = False, add_batch = False, add_pert = True, adata = holdout_adata)
  trained_preds.append(train_pred['latent_x_pred'].X.transpose())
  holdout_preds.append(holdout_pred['latent_x_pred'].X.transpose())

# Make out_dir
Path(out_dir).mkdir(exist_ok = False)
# Save embeddings
np.savetxt(out_dir + "/train-embedding.tsv", latent_outputs["latent_basal"].X.transpose(), delimiter = "\t")
np.savetxt(out_dir + "/holdout-embedding.tsv", basal_predict['latent_z_basal'].X.transpose(), delimiter = "\t")
for idx in np.arange(len(conditions)):
  np.savetxt(out_dir + f"/train-prediction_{conditions[idx]}.tsv", trained_preds[idx], delimiter = "\t")
  np.savetxt(out_dir + f"/holdout-prediction_{conditions[idx]}.tsv", holdout_preds[idx], delimiter = "\t")



session_info.show()
print("Python done")
