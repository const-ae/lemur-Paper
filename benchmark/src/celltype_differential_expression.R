library(SingleCellExperiment)
source("src/utils/config_helper.R")

pa <- argparser::arg_parser("Run differential expression per celltype")
pa <- argparser::add_argument(pa, "--data_id", type = "character", help = "The id of a file in output/results") 
pa <- argparser::add_argument(pa, "--dataset_config", type = "character", nargs = 1, help = "The name of a dataset in datasets.yaml") 
pa <- argparser::add_argument(pa, "--dataset_config_override", type = "character", default = "", nargs = Inf, help = "Override settings from datasets.yaml") 

pa <- argparser::add_argument(pa, "--test_method", type = "character", default = "glmGamPoi", help = "Select the test method")
pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The seed")


pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --data_id 1b927b8a595c6-75ad83f3bbbb6
#                             --dataset_config reyfman
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

config <- get_data_config(pa$dataset_config, pa$dataset_config_override)
out_dir <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

sce <- qs::qread(file.path(pa$working_dir, "results", pa$data_id))

group_var <- colData(sce)[[config$cell_type_column]]
psce_celltype <- glmGamPoi::pseudobulk(sce, group_by = glmGamPoi::vars(fake_condition, group = group_var, sample),
                                       aggregation_functions = list(mod_counts = "rowSums2", .default = "rowMeans2"))

de_res <- purrr::list_rbind(lapply(unique(group_var), \(lvl){
  message("Testing lvl ", lvl)
  if(is.na(lvl)){
    return(NULL)
  }
  de <- tryCatch({
    psce_data <- psce_celltype[,!is.na(psce_celltype$group) & psce_celltype$group == lvl]
    if(pa$test_method == "glmGamPoi"){
      fit <- glmGamPoi::glm_gp(psce_data, ~ fake_condition, size_factors = "ratio", verbose = TRUE)
      glmGamPoi::test_de(fit, contrast = cond(fake_condition = "fake_ctrl") - cond(fake_condition = "fake_trt"))
    }else if(pa$test_method == "edgeR"){
      edger_y <- edgeR::DGEList(counts = assay(psce_data, "counts"))
      edger_design_matrix <- model.matrix(~ fake_condition, colData(psce_data))
      edger_y <- edgeR::calcNormFactors(edger_y)
      edger_y <- edgeR::estimateDisp(edger_y, edger_design_matrix)
      edger_fit <- edgeR::glmQLFit(edger_y, edger_design_matrix, abundance.trend = TRUE, robust = TRUE)
      edger_fit <- edgeR::glmQLFTest(edger_fit, contrast = c(0, -1))
      edger_res <- edgeR::topTags(edger_fit, n = nrow(edger_y), sort.by = "none")$table
      data.frame(name = rownames(edger_res), pval = edger_res$PValue, adj_pval = edger_res$FDR,
                        f_statistic = edger_res$F, df1 = edger_fit$df.test, df2 = edger_fit$df.total, lfc = edger_res$logFC)
    }
  }, error = function(e){
    data.frame(name = rownames(psce_celltype), pval = NA_real_, adj_pval = NA_real_,
               f_statistic = NA_real_, df1 = NA_real_, df2 = NA_real_, lfc = NA_real_)
  })
  de$cell_type <- lvl
  tibble::remove_rownames(de)
}))

# Save everything
dir.create(out_dir)
saveRDS(de_res, file.path(out_dir, "de_results.RDS"))
saveRDS(group_var, file.path(out_dir, "cluster_assignment.RDS"))

#### Session Info
sessionInfo()


