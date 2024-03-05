
suppressWarnings({suppressPackageStartupMessages(library(argparse))})
suppressWarnings({suppressPackageStartupMessages(library(RANN))})
suppressWarnings({suppressPackageStartupMessages(library(data.table))})
suppressWarnings({suppressPackageStartupMessages(library(Seurat))})
suppressWarnings({suppressPackageStartupMessages(library(dplyr))})

loadData <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("File not found.")
  }

  # load data
  data <- data.table::fread(filepath, header = TRUE, sep = "\t",  check.names = FALSE, data.table = FALSE)

  
  row.names(data) <- data[[1]]
  row.names(data) <- gsub(".", "-", row.names(data), fixed = T)
  data[[1]] <- NULL
  
  return(data)
}



smoothDatakNN <- function(output_dir, suffix, max_pcs, seed){

  #message('    Started smoothing binned predicted values by kNN')

  predicted_df <- loadData(paste0(output_dir,'/binned_df',suffix,'.txt'))
  ranked_df <- loadData(paste0(output_dir,'/ranked_df',suffix,'.txt'))
  top_genes <- scan(paste0(output_dir,'/top_var_genes',suffix,'.txt'), what=character(), quiet=TRUE)
  set.seed(seed)

  smoothing_by_KNN <- function(score, umap_coordinates) {
    maxk <- min(30,max(3,round(.005*length(score),0)))
    nearest <- nn2(umap_coordinates, umap_coordinates,k=maxk)
    scaffoldb <- score
    scaffold2 <- sapply(1:length(scaffoldb), function(i) median(scaffoldb[nearest$nn.idx[i,c(1:maxk)]]))
    scaffold2
  }

  suppressMessages({
  suppressWarnings({
    seurat_obj <- CreateSeuratObject(counts = as.matrix(ranked_df), data = as.matrix(ranked_df), min.cells = 0, min.features = 0)
  })
  VariableFeatures(seurat_obj) <- top_genes
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = min(ncol(seurat_obj) - 1, max_pcs))
 
  threshold_var_explained <- 0.5
  stdev_explained <- Stdev(object = seurat_obj, reduction = 'pca')
  var_explained <- stdev_explained * stdev_explained
  var_explained <- data.frame(var_explained)

  var_explained[,'var_explained'] <- var_explained[,'var_explained'] / seurat_obj@reductions$pca@misc$total.variance

  if (sum(var_explained[,'var_explained']) < threshold_var_explained) {
    num_pcs <- min(max_pcs, ncol(seurat_obj) - 1)
  } else {
    num_pcs <- which(cumsum(var_explained[,'var_explained']) > threshold_var_explained)[1]
  }
  num_pcs <- max(2, num_pcs)

  umap_coordinates <- seurat_obj@reductions$pca@cell.embeddings[,1:num_pcs]
  })

  knn_score <- smoothing_by_KNN(predicted_df$preKNN_CytoTRACE2_Score,umap_coordinates[, colnames(umap_coordinates)])
  predicted_df$CytoTRACE2_Score <- knn_score

  # re-map potency categories
  ranges <- seq(0, 1, length.out = 7)
  labels <-
    c(
      'Differentiated',
      'Unipotent',
      'Oligopotent',
      'Multipotent',
      'Pluripotent',
      'Totipotent'
    )

  order_vector <-
    cut(predicted_df$CytoTRACE2_Score,
        breaks = ranges,
        labels = labels)

  predicted_df$CytoTRACE2_Potency <- order_vector


  write.table(predicted_df, paste0(output_dir,"/smoothbykNNresult",suffix,".txt"),sep='\t',quote=F)     
            

}

parser <- ArgumentParser()

parser$add_argument("--output-dir", type="character", default="cytotrace2_results", help="Output directory containing intermediate files")
parser$add_argument("--suffix", type="character", default="", help="Suffix of intermediate files")
parser$add_argument("--max-pcs", type="integer", default=200, help="Integer, indicating the maximum number of principal components to use in the smoothing by kNN step (default is 200)")
parser$add_argument("--seed", type="integer", default=14, help="Integer, specifying the seed for reproducibility in random processes (default is 14).")

args <- parser$parse_args()
smoothDatakNN(output_dir = args$output_dir, suffix = args$suffix, max_pcs = args$max_pcs, seed = args$seed)
                        
