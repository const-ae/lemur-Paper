import os
import sys
print(os.getcwd())
sys.path.insert(0, os.getcwd())  


import argparse
import src.utils.config_helper as ch
from pathlib import Path

import numpy as np
from pandas.api.types import is_numeric_dtype
from scipy import sparse
import scanpy as sc
import scvi
import session_info

parser = argparse.ArgumentParser(description='Run scVI')
parser.add_argument('--data_id', dest='data_id', action='store', required = True, help='The id of a file in output/results')
parser.add_argument("--dataset_config", dest = "dataset_config", required = True, help = "The name of a dataset in datasets.yaml")

parser.add_argument("--working_dir", dest = "working_dir", action='store', required = True, help = "The directory that contains the params, results, scripts etc.")
parser.add_argument("--result_id", dest = "result_id", action='store', required = True, help = "The result_id")
args = parser.parse_args()
# args = parser.parse_args(["--data_id", "337dc3137fe9d-e5f290978d2e2",
#     "--dataset_config", "kang", "--working_dir",
#     "/scratch/ahlmanne/lemur_benchmark", "--result_id", "0"])
print(args)

config = ch.get_data_config(args.dataset_config, '')
out_dir = args.working_dir + "/results/" + args.result_id

# --------------------------------------------------------

adata = sc.read_h5ad(args.working_dir + "/results/" +  args.data_id + "/train.h5ad")
holdout_adata = sc.read_h5ad(args.working_dir + "/results/" +  args.data_id + "/holdout.h5ad")

adata.layers[config['assay_counts']] = sparse.csr_matrix(adata.layers[config['assay_counts']])
holdout_adata.layers[config['assay_counts']] = sparse.csr_matrix(holdout_adata.layers[config['assay_counts']])


# Recommended settings according to harmonization tutorial
model_params = {"n_layers": 2, "n_latent": 30, "gene_likelihood": "nb"}

continuous_batch_cov = is_numeric_dtype(adata.obs[config['main_covariate']])

def prepare_adata(adata):
  if continuous_batch_cov:
    adata.obs[config['main_covariate']] = [str(x) for x in adata.obs[config['main_covariate']]]
    
  scvi.model.SCVI.setup_anndata(
      adata,
      layer=config['assay_counts'],
      batch_key=config['main_covariate']
  )
  return adata
  


prepare_adata(adata)
model = scvi.model.SCVI(adata, **model_params)
if args.dataset_config == "goldfarbmuren":
  # There appears to be a bug when the data size accidentally 
  # leads to a size of 1 for the last batch. Fix it using a
  # different batch_size for Goldfarbmuren. (see also https://discourse.scverse.org/t/solo-scvi-train-error-related-to-batch-size/1591)
  model.train(batch_size = 100)
else:
  model.train()

train_latent_outputs = model.get_latent_representation()
prepare_adata(holdout_adata)
holdout_latent_outputs = model.get_latent_representation(holdout_adata)

def shifted_log_transform(counts, overdispersion = 0.05, pseudo_count = None, minimum_overdispersion = 0.001):
    if pseudo_count is None:
        pseudo_count = 1/(4 * overdispersion)
    size_factors = counts.sum(axis=1)
    size_factors = size_factors / np.exp(np.mean(np.log(size_factors)))
    norm_mat = counts / size_factors.reshape((counts.shape[0], 1))
    overdispersion = 1/(4 * pseudo_count)
    res = 1/np.sqrt(overdispersion) * np.log1p(4 * overdispersion * norm_mat)
    return res

conditions = config['contrast']
if continuous_batch_cov:
  conditions = [str(x) for x in conditions]
  

trained_preds = []
holdout_preds = []
for cond in conditions:
  train_pred  = model.get_normalized_expression(transform_batch = cond, library_size = "latent")
  holdout_pred = model.get_normalized_expression(adata = holdout_adata, transform_batch = cond, library_size = "latent")
   # Apply the same log transformation as to the original counts
  train_pred_mat = shifted_log_transform(train_pred.to_numpy())
  holdout_pred_mat = shifted_log_transform(holdout_pred.to_numpy())
  trained_preds.append(train_pred_mat.transpose())
  holdout_preds.append(holdout_pred_mat.transpose())
  
  


# Make out_dir
Path(out_dir).mkdir(exist_ok = False)
# Save embeddings
np.savetxt(out_dir + "/train-embedding.tsv", train_latent_outputs.transpose(), delimiter = "\t")
np.savetxt(out_dir + "/holdout-embedding.tsv", holdout_latent_outputs.transpose(), delimiter = "\t")
for idx in np.arange(len(conditions)):
  np.savetxt(out_dir + f"/train-prediction_{conditions[idx]}.tsv", trained_preds[idx], delimiter = "\t")
  np.savetxt(out_dir + f"/holdout-prediction_{conditions[idx]}.tsv", holdout_preds[idx], delimiter = "\t")



session_info.show()
print("Python done")
