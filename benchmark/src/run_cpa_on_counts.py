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
# The setup follows the tutorial at https://cpa-tools.readthedocs.io/en/latest/tutorials/Kang.html

parser = argparse.ArgumentParser(description='Run CPA with Kangs parameters')
parser.add_argument('--data_id', dest='data_id', action='store', required = True, help='The id of a file in output/results')
parser.add_argument("--dataset_config", dest = "dataset_config", required = True, help = "The name of a dataset in datasets.yaml")

parser.add_argument("--working_dir", dest = "working_dir", action='store', required = True, help = "The directory that contains the params, results, scripts etc.")
parser.add_argument("--result_id", dest = "result_id", action='store', required = True, help = "The result_id")
args = parser.parse_args()
# args = parser.parse_args(["--data_id", "f1809110732c2-e5f290978d2e2",
#     "--dataset_config", "kang", "--working_dir",
#     "/scratch/ahlmanne/lemur_benchmark", "--result_id", "0"])
print(args)

config = ch.get_data_config(args.dataset_config, '')
out_dir = args.working_dir + "/results/" + args.result_id

# --------------------------------------------------------

adata = sc.read_h5ad(args.working_dir + "/results/" +  args.data_id + "/train.h5ad")
holdout_adata = sc.read_h5ad(args.working_dir + "/results/" +  args.data_id + "/holdout.h5ad")

adata.X = adata.layers[config['assay_counts']].copy()
holdout_adata.X = holdout_adata.layers[config['assay_counts']].copy()

if is_numeric_dtype(adata.obs[config['main_covariate']]):
  adata.obs['key'] = np.where(adata.obs[config['main_covariate']] == config['contrast'][0], "ctrl", "pert")
  holdout_adata.obs['key'] = np.where(holdout_adata.obs[config['main_covariate']] == config['contrast'][0], "ctrl", "pert")
  cpa.CPA.setup_anndata(adata, 
                        perturbation_key='key',
                        control_group = config['contrast'][0], 
                        dosage_key = config['main_covariate'],
                        categorical_covariate_keys=[],
                        is_count_data=True,
                        max_comb_len=1,
                       )
  cpa.CPA.setup_anndata(holdout_adata, 
                        perturbation_key='key',
                        control_group = config['contrast'][0], 
                        dosage_key = config['main_covariate'],
                        categorical_covariate_keys=[],
                        is_count_data=True,
                        max_comb_len=1,
                       )
else:
  adata.obs['dose'] = adata.obs[config['main_covariate']].apply(lambda x: '+'.join(['1.0' for _ in x.split('+')]))
  holdout_adata.obs['dose'] = holdout_adata.obs[config['main_covariate']].apply(lambda x: '+'.join(['1.0' for _ in x.split('+')]))
  cpa.CPA.setup_anndata(adata,
                        perturbation_key=config['main_covariate'],
                        control_group=config['contrast'][-1],
                        dosage_key='dose',
                        categorical_covariate_keys=[],
                        is_count_data=True,
                        max_comb_len=1,
                       )
  cpa.CPA.setup_anndata(holdout_adata,
                        perturbation_key=config['main_covariate'],
                        control_group=config['contrast'][-1],
                        dosage_key='dose',
                        categorical_covariate_keys=[],
                        is_count_data=True,
                        max_comb_len=1,
                       )






model_params = {
    "n_latent": 64,
    "recon_loss": "nb",
    "doser_type": "linear",
    "n_hidden_encoder": 128,
    "n_layers_encoder": 2,
    "n_hidden_decoder": 512,
    "n_layers_decoder": 2,
    "use_batch_norm_encoder": True,
    "use_layer_norm_encoder": False,
    "use_batch_norm_decoder": False,
    "use_layer_norm_decoder": True,
    "dropout_rate_encoder": 0.0,
    "dropout_rate_decoder": 0.1,
    "variational": False,
    "seed": 6977,
}

trainer_params = {
    "n_epochs_kl_warmup": None,
    "n_epochs_pretrain_ae": 30,
    "n_epochs_adv_warmup": 50,
    "n_epochs_mixup_warmup": 0,
    "mixup_alpha": 0.0,
    "adv_steps": None,
    "n_hidden_adv": 64,
    "n_layers_adv": 3,
    "use_batch_norm_adv": True,
    "use_layer_norm_adv": False,
    "dropout_rate_adv": 0.3,
    "reg_adv": 20.0,
    "pen_adv": 5.0,
    "lr": 0.0003,
    "wd": 4e-07,
    "adv_lr": 0.0003,
    "adv_wd": 4e-07,
    "adv_loss": "cce",
    "doser_lr": 0.0003,
    "doser_wd": 4e-07,
    "do_clip_grad": True,
    "gradient_clip_value": 1.0,
    "step_size_lr": 10,
}



model = cpa.CPA(adata, **model_params)
model.train(max_epochs=2000,
            batch_size=512,
            plan_kwargs=trainer_params,
            early_stopping_patience=5,
            check_val_every_n_epoch=5,
           )


# Get integrated embeddings
latent_outputs = model.get_latent_representation(adata, batch_size=1024)
basal_predict = model.custom_predict(basal = True, add_batch = False, add_pert = False, adata = holdout_adata)



def shifted_log_transform(counts, overdispersion = 0.05, pseudo_count = None, minimum_overdispersion = 0.001):
    if pseudo_count is None:
        pseudo_count = 1/(4 * overdispersion)
    size_factors = counts.sum(axis=1)
    size_factors = size_factors / np.exp(np.mean(np.log(size_factors)))
    norm_mat = counts / size_factors.reshape((counts.shape[0], 1))
    overdispersion = 1/(4 * pseudo_count)
    res = 1/np.sqrt(overdispersion) * np.log1p(4 * overdispersion * norm_mat)
    return res

# Make predictions
conditions = config['contrast']
trained_preds = []
holdout_preds = []
for cond in conditions:
  if cond == config['contrast'][0]:
    train_pred = model.custom_predict(basal = True, add_batch = False, add_pert = False)
    holdout_pred = model.custom_predict(basal = True, add_batch = False, add_pert = False, adata = holdout_adata)
  else:
    train_pred = model.custom_predict(basal = False, add_batch = False, add_pert = True)
    holdout_pred = model.custom_predict(basal = False, add_batch = False, add_pert = True, adata = holdout_adata)
    
  # Apply the same log transformation as to the original counts
  train_pred_mat = shifted_log_transform(train_pred['latent_x_pred'].X)
  holdout_pred_mat = shifted_log_transform(holdout_pred['latent_x_pred'].X)
  trained_preds.append(train_pred_mat.transpose())
  holdout_preds.append(holdout_pred_mat.transpose())

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
