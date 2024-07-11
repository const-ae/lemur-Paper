import argparse
from pathlib import Path
import json 

import numpy as np
import scanpy as sc
from gears import PertData, GEARS
import pickle
import session_info


parser = argparse.ArgumentParser(description='Prepare data for combinatorial perturbation prediction')
parser.add_argument('--dataset_name', dest='dataset_name', action='store', required = True, help='The id of a file in output/results')
parser.add_argument('--seed', dest = 'seed', action = 'store', help = "The seed of the run", default = 1, type = int)

parser.add_argument("--working_dir", dest = "working_dir", action='store', required = True, help = "The directory that contains the params, results, scripts etc.")
parser.add_argument("--result_id", dest = "result_id", action='store', required = True, help = "The result_id")
args = parser.parse_args()
# args = parser.parse_args(["--dataset_name", "norman",
#     "--working_dir",
#     "/scratch/ahlmanne/lemur_benchmark", "--result_id", "0"])
print(args)

outfile = args.working_dir + "/results/" + args.result_id
np.random.seed(args.seed)



# -------------------------

pert_data_folder = Path("data/gears_pert_data")

pert_data = PertData(pert_data_folder)
pert_data.load(args.dataset_name)

if args.dataset_name == "norman":
  norman_adata = pert_data.adata
  conds = norman_adata.obs['condition'].values.categories.tolist()
  single_pert = [x for x in conds if 'ctrl' in x]
  double_pert = np.setdiff1d(conds, single_pert).tolist()
  double_training = np.random.choice(double_pert, size=len(double_pert) // 2, replace=False).tolist()
  double_test = np.setdiff1d(double_pert, double_training).tolist()[0]
  double_holdout = np.setdiff1d(double_pert, double_training + [double_test]).tolist()
  
  set2conditions = {
      "train": single_pert + double_training,
      "test": [double_test],
      "val": double_holdout
  }
else:
  # adata = sc.read_h5ad(pert_data_folder / args.dataset_name / "perturb_processed.h5ad")
  # conds = adata.obs['condition'].values.categories.tolist()
  # single_pert = [x for x in conds if 'ctrl' in x]
  # training = np.random.choice(single_pert, size=len(single_pert) // 2, replace=False).tolist()
  # test = np.setdiff1d(single_pert, training).tolist()[0]
  # holdout = np.setdiff1d(single_pert, training + [test]).tolist()
  # set2conditions = {
  #     "train": training,
  #     "test": [test],
  #     "val": holdout
  # }
   pert_data.prepare_split(split = 'simulation', seed = args.seed)
   set2conditions = pert_data.set2conditions


with open(outfile, "w") as outfile: 
    json.dump(set2conditions, outfile)
session_info.show()
print("Python done")
