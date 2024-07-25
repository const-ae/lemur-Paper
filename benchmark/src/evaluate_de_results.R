library(tidyverse)
library(SingleCellExperiment)

pa <- argparser::arg_parser("Evaluate DE results")
pa <- argparser::add_argument(pa, "--data", type = "character", nargs = Inf, help = "The name of the dataset")
pa <- argparser::add_argument(pa, "--vals", type = "character", nargs = Inf, help = "The name of the condition")
pa <- argparser::add_argument(pa, "--method", type = "character", nargs = Inf, help = "The name of the method")
pa <- argparser::add_argument(pa, "--data_seeds", type = "character", nargs = Inf, help = "The name of the method")
pa <- argparser::add_argument(pa, "--de_results_id", type = "character", nargs = Inf, help = "The id of the de_results")
pa <- argparser::add_argument(pa, "--data_results_id", type = "character", nargs = Inf, help = "The id of the data job")

pa <- argparser::add_argument(pa, "--seed", type = "integer", default = 1, nargs = 1, help = "The seed")
pa <- argparser::add_argument(pa, "--nominal_fdr_min", type = "numeric", default = 0.01, nargs = 1, help = "The seed")
pa <- argparser::add_argument(pa, "--nominal_fdr_max", type = "numeric", default = 0.2, nargs = 1, help = "The seed")
pa <- argparser::add_argument(pa, "--large_de_fraction", type = "numeric", default = 2/3, nargs = 1, help = "The seed")
pa <- argparser::add_argument(pa, "--small_de_n", type = "numeric", default = 500, nargs = 1, help = "The seed")
pa <- argparser::add_argument(pa, "--pos_min_cells", type = "numeric", default = 10, nargs = 1, help = "The seed")
pa <- argparser::add_argument(pa, "--pos_fraction", type = "numeric", default = 0.6, nargs = 1, help = "The seed")
pa <- argparser::add_argument(pa, "--neg_fraction", type = "numeric", default = 0.1, nargs = 1, help = "The seed")


pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)
# pa <- argparser::parse_args(pa, argv = r"(
#                             --data angelidis angelidis angelidis
#                             --vals 24m 24m 24m
#                             --method lemur celltype cluster_glmGamPoi
#                             --de_results_id 721ec7bc35156-6770f6509f9dc ef21dda52c1c4-6770f6509f9dc 001192d79f1ba-7a1115ed08c6c
#                             --data_results_id 1b927b8a595c6-aabe6d6194ec7 1b927b8a595c6-aabe6d6194ec7 1b927b8a595c6-aabe6d6194ec7
#                             --working_dir /scratch/ahlmanne/lemur_benchmark
# )" |> stringr::str_trim() |> stringr::str_split("\\s+"))

print(pa)
set.seed(pa$seed)

stopifnot(pa$large_de_fraction > pa$medium_de_fraction)

outdir <- file.path(pa$working_dir, "results/", pa$result_id)
# ---------------------------------------

is_modified_group <- function(genes, is_modified, modified_cells, groups, pos_min_cells = 10, pos_fraction = 0.6, neg_fraction = 0.1){
  stopifnot(length(genes) == length(is_modified))
  stopifnot(length(genes) == length(modified_cells))
  stopifnot(is.logical(is_modified))
  stopifnot(is.list(modified_cells))
  mc_length <- lengths(modified_cells)
  stopifnot(length(setdiff(unique(mc_length), 0)) == 1)
  if(any(is.na(groups))){
    if(is.matrix(groups)) stop("NA's are not allowed in matrix groups")
    filter <- ! is.na(groups)
    modified_cells <- lapply(modified_cells, \(x) if(!is.null(x)){
      x[filter]
    }else NULL)
    groups <- groups[filter]
  }
  
  if(! is.matrix(groups) && ! is(groups, "Matrix")){
    groups <- t(lemur:::one_hot_encoding(groups))
  }
  group_tab <- MatrixGenerics::colSums2(groups, useNames = TRUE)
  group_names <- names(group_tab) %||% seq_along(group_tab)
  
  bind_rows(lapply(seq_along(genes), \(idx){
    gene <- genes[idx]
    mod <- if(is_modified[idx]){
      de_per_group <- MatrixGenerics::colSums2(groups, rows = modified_cells[[idx]], useNames = TRUE)
      group_fraction_de <- de_per_group / as.vector(group_tab)
      case_when(
        rank(-group_fraction_de) == 1 ~ TRUE,   # The highest scoring cluster is always DE
        de_per_group > pos_min_cells & group_fraction_de > pos_fraction ~ TRUE,  # If enough (abs. and rel.) cells in the cluster are DE the cluster is consider DE
        de_per_group < neg_fraction ~ FALSE, # If too few (abs.) cells in the cluster are DE the cluster is consider not DE
        TRUE ~ NA # Otherwise not making a decision
      )
    }else{
      rep(FALSE, length(group_tab))
    } 
    tibble(gene = gene, group = group_names, is_modified = mod)
  }))
}


res <- tibble(data = pa$data,
              vals = pa$vals,
              seeds = pa$data_seeds,
              method = pa$method, 
              de_results_id = pa$de_results_id,
              data_results_id = pa$data_results_id)

de_results <- res %>%
  mutate(de = map2(de_results_id, seq_len(n()), \(id, index){
    message("Read ", id, " at index ", index, "/", nrow(res))
    res <- as_tibble(readRDS(file.path(pa$working_dir, "results", id, "de_results.RDS"))) 
    if("cell_type" %in% colnames(res)){
      arrange(res, name, cell_type)
    }else{
      arrange(res, name)
    }
  }))

data_info <- res %>%
  distinct(data, vals, seeds, data_results_id) %>%
  mutate(n_cells = map_dbl(data_results_id, \(id){
    sce <- qs::qread(file.path(pa$working_dir, "results", id))
    ncol(sce)
  }))

gene_info <- res %>%
  distinct(data, vals, seeds, data_results_id) %>%
  mutate(genes = map(data_results_id, \(id){
    as_tibble(SummarizedExperiment::rowData(qs::qread(file.path(pa$working_dir, "results", id)))) %>%
      mutate(de_size = map_dbl(is_de_cell, \(x) sum(x != 0)))
  })) %>%
  dplyr::select(-data_results_id)

global_ground_truths <- res %>%
  filter(str_starts(method, "lemur") | str_starts(method, "global")) %>%
  left_join(gene_info, by = c("data", "vals", "seeds")) %>%
  mutate(ground_truth = map(genes, \(x) x |> transmute(gene = name, group = NA, is_modified = is_simulated) |> arrange(gene)))

cluster_ground_truths <- res %>%
  filter(str_starts(method, "cluster") | str_starts(method, "celltype") | str_starts(method, "miloDE")) %>%
  left_join(gene_info, by = c("data", "vals", "seeds")) %>%
  mutate(ground_truth = map2(genes, de_results_id, \(x, y){
    is_modified_group(x$name, is_modified = x$is_simulated, x$is_de_cell, readRDS(file.path(pa$working_dir, "results", y, "cluster_assignment.RDS")),
                      pos_min_cells = pa$pos_min_cells, pos_fraction = pa$pos_fraction, neg_fraction = pa$neg_fraction) %>%
      arrange(gene, group)
  }))  

total_de_info <- left_join(
  de_results,
  bind_rows(global_ground_truths, cluster_ground_truths) %>% dplyr::select(-c(de_results_id, data_results_id)), 
  by = c("data", "vals", "seeds", "method")) %>%
  left_join(data_info, by = c("data", "vals", "seeds", 'data_results_id'))

nominal_fdr <- c(1e-3, seq(pa$nominal_fdr_min, pa$nominal_fdr_max, by = 0.01))
power_res <- total_de_info  %>%
  summarize(results = map2(de, ground_truth, \(de, gt){
    merged <- if("cell_type" %in% colnames(de)){
      inner_join(de %>% mutate(cell_type = as.character(cell_type)), gt, by = c("name" = "gene", "cell_type" = "group"))
    }else{
      inner_join(de, gt, by = c("name" = "gene"))
    }
    bind_rows(lapply(nominal_fdr, \(thres){
      merged %>%
        mutate(is_signif = ! is.na(adj_pval) & adj_pval < thres) %>%
        summarize(TPR = sum(is_signif & is_modified, na.rm = TRUE) / sum(is_modified, na.rm = TRUE),
                  TP = sum(is_signif & is_modified, na.rm = TRUE),
                  FDR = sum(is_signif & !is_modified, na.rm = TRUE) / sum(is_signif, na.rm = TRUE)) %>%
        mutate(nominal_fdr = thres)
    }))
  }),
  strat_results =  pmap(list(de, ground_truth, genes, n_cells), \(de, gt, ge, n_cells){
    merged <- (if("cell_type" %in% colnames(de)){
      inner_join(de %>% mutate(cell_type = as.character(cell_type)), gt, by = c("name" = "gene", "cell_type" = "group"))
    }else{
      inner_join(de, gt, by = c("name" = "gene"))
    }) %>%
      inner_join(ge %>% transmute(name, de_size), by = c("name")) %>%
      filter(de_size > 0) %>%
      mutate(de_size_cat = case_when(
        de_size > pa$large_de_fraction * .env$n_cells ~ "large",
        de_size <= pa$small_de_n ~ "small",
        TRUE ~ "medium"
      ))
    
    bind_rows(lapply(nominal_fdr, \(thres){
      merged %>%
        mutate(is_signif = ! is.na(adj_pval) & adj_pval < thres) %>%
        summarize(TPR = sum(is_signif & is_modified, na.rm = TRUE) / sum(is_modified, na.rm = TRUE),
                  TP = sum(is_signif & is_modified, na.rm = TRUE),
                  .by = de_size_cat) %>%
        mutate(nominal_fdr = thres)
    }))
  }),
  .by = c(data, vals, seeds, method)) 


# Save everything
tmp_out_dir <- paste0(outdir, "-tmp")
dir.create(tmp_out_dir)
# This takes hours to write
# saveRDS(total_de_info, file.path(tmp_out_dir, "de_results.RDS"))
saveRDS(power_res, file.path(tmp_out_dir, "fdr_tpr_results.RDS"))
file.rename(tmp_out_dir, outdir)

#### Session Info
sessionInfo()


