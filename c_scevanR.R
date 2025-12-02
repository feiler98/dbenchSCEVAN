library(Seurat)
library(SCEVAN)
library(glue)

csv_to_seurat <- function(path_rcm) {
  raw_rcm <- read.csv(path_rcm, header = T)
  rownames(raw_rcm) <-  make.names(raw_rcm$Gene, unique = TRUE)
  raw_rcm$Gene <- NULL
  seurat_obj <- CreateSeuratObject(counts = raw_rcm, assay="RNA")
  return (GetAssayData(seurat_obj ,assay = 'RNA', layer = 'counts'))
}

# norm_cells is a vector here!
r_run_scevan <- function(path_file,
                         sample_tag,
                         par_cores=10,
                         norm_cell=NULL,
                         ngenes_chr=5,
                         perc_genes=10,
                         beta_vega=0.5) {
    raw_seurat <- csv_to_seurat(path_file)
    results <- pipelineCNA(raw_seurat,
                           sample=sample_tag,
                           par_cores=par_cores,
                           norm_cell=norm_cell,
                           ngenes_chr=ngenes_chr,
                           perc_genes=perc_genes,
                           beta_vega=beta_vega)
    write.csv(results, glue("{sample_tag}__results.csv"))
}