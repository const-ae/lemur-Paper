library(SingleCellExperiment)
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Run differential using miloDE")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--dataset_config", type = "character", nargs = 1, help = "The name of a dataset in datasets.yaml") 
pa <- argparser::add_argument(pa, "--dataset_config_override", type = "character", default = "", nargs = Inf, help = "Override settings from datasets.yaml") 

pa <- argparser::add_argument(pa, "--padjust_method", type = "character", default = "across_genes", nargs = 1, help = "Which method to use for the p-value adjustment")
pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The seed")


pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --data_id 060e8c8f52b5c-5f4f9c6b45716
#                             --dataset_config angelidis
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

adjustment_method <- match.arg(pa$padjust_method, c("across_genes", "across_nhoods", "across_all"))

config <- get_data_config(pa$dataset_config, pa$dataset_config_override)
out_dir <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

sce <- qs::qread(file.path(pa$working_dir, "results", pa$data_id))

hvg <- order(-rowVars(assay(sce, config$assay_continuous)))
pca <- lemur:::pca(assay(sce, config$assay_continuous)[hvg[seq_len(min(nrow(sce), 500))],], n = min(nrow(sce), ncol(sce), 50))
harm <- harmony::RunHarmony(pca$embedding, vars_use = c("fake_condition", "sample"), meta_data = colData(sce), lambda = c(1,1))

milo_sce <- sce
assay(milo_sce, "counts") <- assay(milo_sce, config$assay_counts)
assay(milo_sce, "logcounts") <- assay(milo_sce, config$assay_continuous)
reducedDim(milo_sce, "harmony") <- harm
milo_sce <- miloDE::assign_neighbourhoods(milo_sce, "harmony")
system.time(
  milo_de <- miloDE::de_test_neighbourhoods(milo_sce, sample_id = "sample", design = ~ fake_condition, covariates = "fake_condition")
)
milo_de$adj_pval <- if(adjustment_method == "across_genes"){
  milo_de$pval_corrected_across_genes
}else if(adjustment_method == "across_nhoods"){
  milo_de$pval_corrected_across_nhoods
}else if(adjustment_method == "across_all"){
  p.adjust(milo_de$pval, "BH")
}


de_res <- data.frame(name = milo_de$gene, pval = milo_de$pval, adj_pval = milo_de$adj_pval,
                     lfc = milo_de$logFC, cell_type = as.character(milo_de$Nhood_center))

# Save everything
dir.create(out_dir)
saveRDS(de_res, file.path(out_dir, "de_results.RDS"))
saveRDS(miloR::nhoods(milo_sce), file.path(out_dir, "cluster_assignment.RDS"))

#### Session Info
sessionInfo()


