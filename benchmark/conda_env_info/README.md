This directory contains the details of the packages in each of the conda 
environments used for the benchmark.

The files were produced by calling
```
conda env export -n pylemur_env > conda_env_info/pylemur_env.yml
```

You can recreate the conda environment by calling
```
conda env create -f path/to/environment.yml
```