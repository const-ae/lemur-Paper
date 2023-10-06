

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
  
  
  list(angelidis = get_angelidis_data, kang = get_kang_data, skinnider = get_skinnider_data)
})()
