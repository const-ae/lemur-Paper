#!/bin/bash

# Remove old folder content
rm -r benchmark
# Make new folder
mkdir benchmark
mkdir benchmark/renv
mkdir benchmark/data

# Copy relevant files (i.e., exclude renv library because it is too big)
scp -r ahlmanne@seneca:~/projects/lemur-Paper-benchmark/{README.md,renv.lock,.Rprofile,run_benchmarks.sh,conda_env_package_versions.txt} benchmark/.
scp -r ahlmanne@seneca:~/projects/lemur-Paper-benchmark/renv/activate.R benchmark/renv/activate.R
scp -r ahlmanne@seneca:~/projects/lemur-Paper-benchmark/submission benchmark/submission
scp -r ahlmanne@seneca:~/projects/lemur-Paper-benchmark/src benchmark/src
scp -r ahlmanne@seneca:~/projects/lemur-Paper-benchmark/output benchmark/output
scp -r ahlmanne@seneca:~/projects/lemur-Paper-benchmark/data/{de_data_spec.tsv,datasets.yaml,.gitignore} benchmark/data/.