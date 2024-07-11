from pathlib import Path
import numpy as np
import pandas as pd
import json
import tempfile
import shutil
import argparse
import session_info

import scanpy as sc
import pylemur

parser = argparse.ArgumentParser(description='Run pylemur additive model')
parser.add_argument('--dataset_name', dest='dataset_name', action='store', required = True, help='The id of a file in output/results')
parser.add_argument('--test_train_config_id', dest = 'test_train_config_id', action = 'store', required = True, help = "The ID of the test/train/holdout run")
parser.add_argument('--align_harmony', dest = 'align_harmony', action = 'store', help = "Should it run LEMUR's align_with_harmony function", default = False, type = bool)
parser.add_argument('--n_embedding', dest = 'n_embedding', action = 'store', help = "The number of latent dimensions for LEMUR", default = 2, type = int)

parser.add_argument('--seed', dest = 'seed', action = 'store', help = "The seed of the run", default = 1, type = int)

parser.add_argument("--working_dir", dest = "working_dir", action='store', required = True, help = "The directory that contains the params, results, scripts etc.")
parser.add_argument("--result_id", dest = "result_id", action='store', required = True, help = "The result_id")

args = parser.parse_args()
# args = parser.parse_args(["--dataset_name", "norman",
#     "--test_train_config_id", "8443ed21d2ac4-f8716281f960b", "--working_dir",
#     "/scratch/ahlmanne/lemur_benchmark", "--result_id", "0"])
print(args)

out_dir = args.working_dir + "/results/" + args.result_id

np.random.seed(args.seed)
# --------------------------------------------------------

pert_data_folder = Path("data/gears_pert_data/")
adata = sc.read_h5ad(pert_data_folder / args.dataset_name / "perturb_processed.h5ad")
adata.layers["logcounts"] = adata.X.copy()

with open(args.working_dir + "/results/" + args.test_train_config_id) as json_file:
  set2conditions = json.load(json_file)

all_valid_conds = set(sum(set2conditions.values(), []))
adata = adata[[x in all_valid_conds for x in adata.obs.condition],:]

def comb_to_design_matrix(lst, levels = None):
  if levels is None:
    levels = set()
    for x in lst:
      levels.update(x)
    levels = list(levels)
    levels.sort()
  if not isinstance(levels, np.ndarray):
    levels = np.array(levels)
  if len(levels.shape) != 1:
    raise ValueError("Levels must be a one dimensional array")
  
  nrow = len(lst)
  design_mat = np.zeros((nrow, len(levels)))
  for cell_idx in range(nrow):
    for cond in lst[cell_idx]:
      col_idx = np.where(levels == cond)    
      design_mat[cell_idx, col_idx] = 1
  return pd.DataFrame(design_mat, columns=levels)

training_df = (pd.DataFrame({"training": set2conditions.keys(),
              "condition": set2conditions.values()})
  .explode("condition"))

adata.obs = adata.obs.merge(training_df, how = "left", on = "condition")

all_split_conds = [x.split("+") for x in adata.obs["condition"]]
design_mat = comb_to_design_matrix(all_split_conds)
design_mat['ctrl'] = 1

training_subset = adata.obs['training'] == "train"
model = pylemur.tl.LEMUR(adata[training_subset,:], design = design_mat[training_subset], n_embedding= args.n_embedding)
model.fit()
if args.align_harmony:
  model.align_with_harmony()

conds = list(set(adata.obs["condition"].tolist()))
split_conds = [x.split("+") for x in conds]
center = model.embedding.mean(axis = 0).reshape(1, -1)
preds = [model.predict(embedding = center, new_condition=comb_to_design_matrix([c + ["ctrl"]], levels=design_mat.columns))[0,:] for c in split_conds]

all_pred_vals = dict(zip(conds, preds))
all_pred_vals = {k: v.tolist() for k, v in all_pred_vals.items()}

tmp_out_dir = tempfile.mkdtemp()
with open(f"{tmp_out_dir}/all_predictions.json", 'w', encoding="utf8") as handle:
    json.dump(all_pred_vals, handle, indent = 4)
with open(f"{tmp_out_dir}/gene_names.json", 'w', encoding="utf8") as handle:
    json.dump(adata.var["gene_name"].values.tolist(), handle, indent = 4)

# Move results to out_dir
shutil.move(tmp_out_dir, out_dir)



session_info.show()
print("Python done")
