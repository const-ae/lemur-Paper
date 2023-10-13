

data_loader <- (function(){
  
  library(SingleCellExperiment)
  
  # Angelidis
  get_angelidis_data <- function(){
    if(file.exists("data/angelidis.qs")){
      qs::qread("data/angelidis.qs")
    }else{
      library(tidyverse)
      library(MatrixGenerics)
      if(! file.exists("data/angelidis/mouse_lung_single_cell.RData")){
        dir.create("data/angelidis")
        download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE124872&format=file&file=GSE124872%5Fraw%5Fcounts%5Fsingle%5Fcell%2ERData%2Egz", "data/angelidis/mouse_lung_single_cell.RData.gz")
        R.utils::gunzip("data/angelidis/mouse_lung_single_cell.RData.gz")
      }
      if(! file.exists("data/angelidis/mouse_lung_single_cell_metadata.csv")){
        download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE124872&format=file&file=GSE124872%5FAngelidis%5F2018%5Fmetadata%2Ecsv%2Egz", "data/angelidis/mouse_lung_single_cell_metadata.csv.gz")
        R.utils::gunzip("data/angelidis/mouse_lung_single_cell_metadata.csv.gz")
      }
      load("data/angelidis/mouse_lung_single_cell.RData") # creates the 'raw_counts' variable
      metadat <- read_csv("data/angelidis/mouse_lung_single_cell_metadata.csv")
      col_data <- metadat %>%
        transmute(barcode = str_split_fixed(...1, ":", 3)[,3],
                  mouse_id = identifier, age = grouping, cell_type = celltype) %>%
        mutate(cell_id = paste0(mouse_id, ":", barcode))
      tmp <- tibble(X1 = colnames(raw_counts)) %>%
        separate(X1, c("mouse_id", "barcode"), sep = ":")
      stopifnot(tmp$barcode == col_data$barcode)
      colnames(raw_counts) <- col_data$cell_id
      colnames(raw_counts) <- col_data$cell_id
      sce <- SingleCellExperiment(list(counts = raw_counts), colData = col_data)
      total_umi <- colSums2(assay(sce, "counts"))
      sce <- sce[rowSums2(assay(sce, "counts")) > 0, total_umi > median(total_umi) / 10 & total_umi < median(total_umi) * 10]
      assay(sce, "logcounts") <- transformGamPoi::shifted_log_transform(sce)
      
      qs::qsave(sce, "data/angelidis.qs")
      sce
    }
  } 
  
  get_aztekin_data <- function(){
    if(file.exists("data/aztekin.qs")){
      qs::qread("data/aztekin.qs")
    }else{
      sce <- scRNAseq::AztekinTailData()
      assay(sce, "logcounts") <- transformGamPoi::shifted_log_transform(sce)
      qs::qsave(sce, "data/aztekin.qs")
      sce
    }
  }
  
  get_bunis_data <- function(){
    if(file.exists("data/bunis.qs")){
      qs::qread("data/bunis.qs")
    }else{
      sce <- scRNAseq::BunisHSPCData()
      assay(sce, "logcounts") <- transformGamPoi::shifted_log_transform(sce)
      colData(sce)$DevStageScoring <- NULL
      qs::qsave(sce, "data/bunis.qs")
      sce
    }
  }
  
  get_goldfarbmuren_data <- function(){
    if(file.exists("data/goldfarbmuren.qs")){
      qs::qread("data/goldfarbmuren.qs")
    }else{
      if(! file.exists("data/goldfarbmuren/GSE134174_Processed_invivo_raw.txt.gz")){
        dir.create("data/goldfarbmuren")
        download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE134174&format=file&file=GSE134174%5FProcessed%5Finvivo%5Fraw%2Etxt%2Egz", "data/goldfarbmuren/GSE134174_Processed_invivo_raw.txt.gz")
      }
      if(! file.exists("data/goldfarbmuren/metadata.tsv.gz")){
        download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE134174&format=file&file=GSE134174%5FProcessed%5Finvivo%5Fmetadata%2Etxt%2Egz", "data/goldfarbmuren/metadata.tsv.gz")
      }
      library(tidyverse)
      library(MatrixGenerics)
      counts <- as(as.matrix(data.table::fread("data/goldfarbmuren/GSE134174_Processed_invivo_raw.txt.gz"), rownames = 1), "dgCMatrix")
      col_data <- readr::read_tsv("data/goldfarbmuren/metadata.tsv.gz")
      stopifnot(all(colnames(counts) == col_data$Cell))
      
      sce <- SingleCellExperiment(list(counts = counts), colData = as.data.frame(col_data))
      sce <- sce[,sce$Smoke_status %in% c("heavy", "never")]
      assay(sce, "logcounts") <- transformGamPoi::shifted_log_transform(sce)
      qs::qsave(sce, "data/goldfarbmuren.qs")
      sce
    }
  }
  
  get_hrvatin_data <- function(){
    if(file.exists("data/hrvatin.qs")){
      qs::qread("data/hrvatin.qs")
    }else{
      if(! file.exists("data/hrvatin/GSE102827_merged_all_raw.csv.gz")){
        dir.create("data/hrvatin")
        options(timeout=1000)
        download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE102827&format=file&file=GSE102827%5Fmerged%5Fall%5Fraw%2Ecsv%2Egz", "data/hrvatin/GSE102827_merged_all_raw.csv.gz")
      }
      if(! file.exists("data/hrvatin/GSE102827_cell_type_assignments.csv.gz")){
        download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE102827&format=file&file=GSE102827%5Fcell%5Ftype%5Fassignments%2Ecsv%2Egz", "data/hrvatin/GSE102827_cell_type_assignments.csv.gz")
      }
      library(tidyverse)
      library(MatrixGenerics)
      counts <- as(as.matrix(data.table::fread("data/hrvatin/GSE102827_merged_all_raw.csv.gz"), rownames = 1), "dgCMatrix")
      col_data <- read_csv("data/hrvatin/GSE102827_cell_type_assignments.csv.gz") %>%
        dplyr::rename(Cell = ...1)
      stopifnot(all(colnames(counts) == col_data$Cell))
      
      sce <- SingleCellExperiment(list(counts = counts), colData = as.data.frame(col_data))
      assay(sce, "logcounts") <- transformGamPoi::shifted_log_transform(sce)
      qs::qsave(sce, "data/hrvatin.qs")
      sce
    }
  }

  get_jakel_data <- function(){
    if(file.exists("data/jakel.qs")){
      qs::qread("data/jakel.qs")
    }else{
      library(tidyverse)
      library(MatrixGenerics)
      if(! file.exists("data/jakel/GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz")){
        dir.create("data/jakel")
        options(timeout=1000)
        download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118257&format=file&file=GSE118257%5FMSCtr%5FsnRNA%5FExpressionMatrix%5FR%2Etxt%2Egz", "data/jakel/GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz")
      }
      if(! file.exists("data/jakel/GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz")){
        download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118257&format=file&file=GSE118257%5FMSCtr%5FsnRNA%5FFinalAnnotationTable%2Etxt%2Egz", "data/jakel/GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz")
      }
      
      counts <- as(as.matrix(data.table::fread("data/jakel/GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz"), rownames = 1), "dgCMatrix")
      col_data <- read_tsv("data/jakel/GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz") %>%
        dplyr::rename(Cell = Detected, detected_genes = genes)
      stopifnot(all(colnames(counts) == col_data$Cell))
      
      sce <- SingleCellExperiment(list(counts = counts), colData = as.data.frame(col_data))
      assay(sce, "logcounts") <- transformGamPoi::shifted_log_transform(sce)
      qs::qsave(sce, "data/jakel.qs")
      sce
    }
  }
  
  get_sathyamurthy_data <- function(){
    if(file.exists("data/sathyamurthy.qs")){
      qs::qread("data/sathyamurthy.qs")
    }else{
      library(tidyverse)
      library(MatrixGenerics)
      if(! file.exists("data/sathyamurthy/GSE103892_Expression_Count_Matrix.txt.gz")){
        dir.create("data/sathyamurthy")
        options(timeout=1000)
        download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103892&format=file&file=GSE103892%5FExpression%5FCount%5FMatrix%2Etxt%2Egz",
                      "data/sathyamurthy/GSE103892_Expression_Count_Matrix.txt.gz")
      }
      if(! file.exists("data/sathyamurthy/GSE103892_Sample_Cell_Cluster_Information.txt.gz")){
        download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103892&format=file&file=GSE103892%5FSample%5FCell%5FCluster%5FInformation%2Etxt%2Egz",
                      "data/sathyamurthy/GSE103892_Sample_Cell_Cluster_Information.txt.gz")
      }
      
      counts <- as(as.matrix(data.table::fread("data/sathyamurthy/GSE103892_Expression_Count_Matrix.txt.gz"), rownames = 1), "dgCMatrix")
      col_data <- read_tsv("data/sathyamurthy/GSE103892_Sample_Cell_Cluster_Information.txt.gz", comment = "TYPE") %>%
        separate(col = "sample_cellbarcode", into = c("sample", "barcode"), remove = FALSE) %>%
        mutate(condition = str_match(sample, "^(\\w+?)\\d{0,2}$")[,2])
      table(col_data$condition)
      # Reorder the rows of col data to match the counts matrix
      col_data <- tibble(sample_cellbarcode = colnames(counts)) %>%
        tidylog::left_join(col_data)
      stopifnot(all(colnames(counts) == col_data$sample_cellbarcode))
      
      sce <- SingleCellExperiment(list(counts = counts), colData = as.data.frame(col_data))
      sce <- sce[,! is.na(sce$cell.type) & sce$cell.type != "discarded" & sce$condition %in% c("form", "rotarod")]
      assay(sce, "logcounts") <- transformGamPoi::shifted_log_transform(sce)
      qs::qsave(sce, "data/sathyamurthy.qs")
      sce
    }
  }
  
  # Download Data from Kang et al.
  get_kang_data <- function(){
    if(file.exists("data/kang.qs")){
      qs::qread("data/kang.qs")
    }else{
      library(tidyverse)
      library(MatrixGenerics)
      if(! file.exists("data/kang_augur_data/augur_protocol_zenodo.tar.gz")){
        options(timeout=1000)
        download.file("https://zenodo.org/record/4473025/files/augur_protocol_zenodo.tar.gz?download=1", "data/kang_augur_data/augur_protocol_zenodo.tar.gz")
      }
      if(! file.exists("data/kang_augur_data/augur_protocol_zenodo")){
        untar("data/kang_augur_data/augur_protocol_zenodo.tar.gz", exdir = "data/kang_augur_data")
      }
      seurat_dat <- readRDS("data/kang_augur_data/augur_protocol_zenodo/rnaseq/processed/Kang2018.rds")
      sce <- Seurat::as.SingleCellExperiment(seurat_dat)
      altExp(sce, "SCT") <- NULL
      altExp(sce, "integrated") <- NULL
      assay(sce, "logcounts") <- transformGamPoi::shifted_log_transform(sce)
      
      qs::qsave(sce, "data/kang.qs")
      sce
    }
  }
  
  get_skinnider_data <- function(){
    if(file.exists("data/skinnider.qs")){
      qs::qread("data/skinnider.qs")
    }else{
      library(tidyverse)
      library(MatrixGenerics)
      if(! file.exists("data/kang_augur_data/augur_protocol_zenodo.tar.gz")){
        options(timeout=1000)
        download.file("https://zenodo.org/record/4473025/files/augur_protocol_zenodo.tar.gz?download=1", "data/kang_augur_data/augur_protocol_zenodo.tar.gz")
      }
      if(! file.exists("data/kang_augur_data/augur_protocol_zenodo")){
        untar("data/kang_augur_data/augur_protocol_zenodo.tar.gz", exdir = "data/kang_augur_data")
      }
      mats <- readRDS("data/kang_augur_data/augur_protocol_zenodo/rnaseq/processed/Skinnider2020_mats.rds")
      meta <- readRDS("data/kang_augur_data/augur_protocol_zenodo/rnaseq/processed/Skinnider2020_meta.rds") %>%
        ungroup()
      sce <- SingleCellExperiment(list(counts = mats$spliced), colData = as.data.frame(meta))
      assay(sce, "logcounts") <- transformGamPoi::shifted_log_transform(sce)
      
      qs::qsave(sce, "data/skinnider.qs")
      sce
    }
  }
  
  get_bhattacherjee_data <- function(){
    # https://doi.org/10.1038/s41467-019-12054-3
    if(file.exists("data/bhattacherjee.qs")){
      qs::qread("data/bhattacherjee.qs")
    }else{
      library(tidyverse)
      library(MatrixGenerics)
      if(! file.exists("data/kang_augur_data/augur_protocol_zenodo.tar.gz")){
        options(timeout=1000)
        download.file("https://zenodo.org/record/4473025/files/augur_protocol_zenodo.tar.gz?download=1", "data/kang_augur_data/augur_protocol_zenodo.tar.gz")
      }
      if(! file.exists("data/kang_augur_data/augur_protocol_zenodo")){
        untar("data/kang_augur_data/augur_protocol_zenodo.tar.gz", exdir = "data/kang_augur_data")
      }
      seurat_dat <- readRDS("data/kang_augur_data/augur_protocol_zenodo/rnaseq/processed/Bhattacherjee2019.rds")
      sce <- Seurat::as.SingleCellExperiment(seurat_dat)
      qs::qsave(sce, "data/bhattacherjee.qs")
      sce
    }
  }
  
  get_squair_suppl_data <- function(){
    if(! file.exists("data/squair_18_single_cell_datasets.tar.gz")){
      options(timeout=2000)
      download.file("https://zenodo.org/record/5048449/files/sc_rnaseq.tar.gz?download=1", "data/squair_18_single_cell_datasets.tar.gz")
    }
    if(! file.exists("data/squair_18_single_cell_datasets")){
      untar("data/squair_18_single_cell_datasets.tar.gz", exdir = "data/squair_18_single_cell_datasets")
    }
    "data/squair_18_single_cell_datasets"
  }
  
  
  get_canogamez_data <- function(){
    # Cano-Gamez, Eddie, Blagoje Soskic, Theodoros I. Roumeliotis, Ernest So, Deborah J. Smyth, Marta Baldrighi, David Willé, et al. “Single-Cell Transcriptomics Identifies an Effectorness Gradient Shaping the Response of CD4+ T Cells to Cytokines.” Nature Communications 11, no. 1 (April 14, 2020): 1801. https://doi.org/10.1038/s41467-020-15543-y.
    # This dataset is cool because there is an effectorness gradient which
    # depends on the underlying cell type (memory or helper T cells).
    # The original publication contains a whole figure about interactions
    if(file.exists("data/canogamez.qs")){
      qs::qread("data/canogamez.qs")
    }else{
      library(tidyverse)
      library(MatrixGenerics)
      dir <- get_squair_suppl_data()
      sce <- do.call(cbind, lapply(list.files("data/squair_18_single_cell_datasets/sc_rnaseq/rds", pattern = "^CanoGamez.*\\.rds", full.names = TRUE), \(file){
        message("Read ", file)
        Seurat::as.SingleCellExperiment(readRDS(file))
      }))
      sce <-  sce[, sce$label != "UNS"]
      assay(sce, "logcounts") <- transformGamPoi::shifted_log_transform(sce)
      qs::qsave(sce, "data/canogamez.qs")
      sce
    }
  }
  
  get_reyfman_data <- function(){
    if(file.exists("data/reyfman.qs")){
      qs::qread("data/reyfman.qs")
    }else{
      library(tidyverse)
      library(MatrixGenerics)
      dir <- get_squair_suppl_data()
      sce <- do.call(cbind, lapply(list.files("data/squair_18_single_cell_datasets/sc_rnaseq/rds", pattern = "^Reyfman.*\\.rds", full.names = TRUE), \(file){
        message("Read ", file)
        Seurat::as.SingleCellExperiment(readRDS(file))
      }))
      assay(sce, "logcounts") <- transformGamPoi::shifted_log_transform(sce)
      qs::qsave(sce, "data/reyfman.qs")
      sce
    }
  }
  
  get_miloDE_mouse_gastrulation_data <- function(){
    if(file.exists("data/mouse_gastrulation.qs")){
      qs::qread("data/mouse_gastrulation.qs")
    }else{
      library(tidyverse)
      library(MatrixGenerics)
      sce_chimera <- MouseGastrulationData::Tal1ChimeraData()
      # select tomato positive cells (tal1 KO); delete doublets and stripped
      sce_chimera = sce_chimera[, !sce_chimera$celltype.mapped %in% c("Doublet" , "Stripped")]
      # delete row for tomato
      sce_chimera = sce_chimera[rownames(sce_chimera) != "tomato-td" , ]
      
      colnames(sce_chimera) <- paste0("chimera_" , colnames(sce_chimera) )
      sce_chimera$cell <- paste0("chimera_" ,sce_chimera$cell)
      sce_chimera$sample <- paste0("chimera_" , sce_chimera$sample )
      sce_chimera$type <- "chimera"
      sce_chimera$tal1 <- ifelse(sce_chimera$tomato , 0 , 1)
      
      samples_e8.5 <- MouseGastrulationData::AtlasSampleMetadata$sample[MouseGastrulationData::AtlasSampleMetadata$stage == "E8.5"]
      
      sce_ref <- do.call(cbind, lapply(samples_e8.5 , function(sample){
        out = EmbryoAtlasData(samples = sample)
        return(out)
      }))
      sce_ref <- sce_ref[, !sce_ref$stripped & !sce_ref$doublet]
      sce_ref$type <- "wt"
      sce_ref$tal1 <- 1
      
      stopifnot(all(rownames(sce_ref) == rownames(sce_chimera)))
      sce <- SingleCellExperiment(list(counts = cbind(counts(sce_ref) , counts(sce_chimera))))
      sce$cell <- c(sce_ref$cell , sce_chimera$cell)
      sce$celltype <- c(sce_ref$celltype , sce_chimera$celltype.mapped)
      sce$sample <- c(sce_ref$sample , sce_chimera$sample)
      sce$type <- c(sce_ref$type , sce_chimera$type)
      sce$tal1 <- c(sce_ref$tal1 , sce_chimera$tal1)
      
      rowData(sce) <- rowData(sce_ref)
      sce <- scuttle::logNormCounts(sce)
      
      qs::qsave(sce, "data/mouse_gastrulation.qs")
      sce
    }
  }
  
  list(angelidis = get_angelidis_data, aztekin = get_aztekin_data, bunis = get_bunis_data,
       goldfarbmuren = get_goldfarbmuren_data,  hrvatin = get_hrvatin_data, jakel = get_jakel_data, 
       sathyamurthy = get_sathyamurthy_data, kang = get_kang_data, skinnider = get_skinnider_data,
       bhattacherjee = get_bhattacherjee_data, canogamez = get_canogamez_data, reyfman = get_reyfman_data,
       mouse_gastrulation = get_miloDE_mouse_gastrulation_data)
})()
