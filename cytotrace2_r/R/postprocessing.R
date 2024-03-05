#' Dispersion Function
#'
#' This function calculates a dispersion metric for a given numeric vector.
#'
#' @param x A numeric vector for which the dispersion metric needs to be calculated.
#' @return A numeric value representing the dispersion metric of the input vector.
#' @export

disp_fn <- function(x) {
  if (length(unique(x)) == 1) {
    return(0)
  }
  else {
    return(var(x)/mean(x))
  }
}

#' Smooth Data with kNN
#'
#' This applies k-nearest neighbor smoothing to datasets with >100 cells.
#' It identifies the nearest neighbors of each cell based on top principal components,
#' and assigns the median score of the cell and its  nearest neighbors, as the final
#' potency score for the cell, and updates potency category accordingly.
#'
#' @param ranked_df dataframe of ranked gene expression (ranked data, transposed) where rows are genes and columns are cells
#' @param predicted_df dataframe of 4 columns containing CytoTRACE 2 model predicted and postprocessed potency scores and categories
#' @param top_genes list of the top 1000 most variable genes (output of dispersion function)
#' @param max_pcs integer indicating the maximum number of principal components to use in the smoothing by kNN step (default is 200).
#' @param seed integer specifying the seed for reproducibility in random processes (default is 14).
#'
#' @return A dataframe with final updated predictions
#' @importFrom Seurat CreateSeuratObject ScaleData RunPCA VariableFeatures Stdev
#' @importFrom RANN nn2
#' @export
#'
#' @examples
#'  ranked_data --> first element in the output of preprocessData()
#'  predicted_df --> output of binData()
#'  ranked_df <- data.frame(t(ranked_data))
#'  rownames(ranked_df) <- row.names(expression_data)
#'  colnames(ranked_df) <- colnames(expression_data)
#'  predicted_df <- smoothDatakNN(ranked_df, predicted_df, top_genes)


smoothDatakNN <- function(ranked_df, predicted_df, top_genes, max_pcs, seed){

  # message('cytotrace2: Started smoothing binned predicted values by kNN')

  set.seed(seed)

  smoothing_by_KNN <- function(score, umap_coordinates) {
    maxk <- min(30,max(3,round(.005*length(score),0)))
    nearest <- RANN::nn2(umap_coordinates, umap_coordinates,k=maxk)
    scaffoldb <- score
    scaffold2 <- sapply(1:length(scaffoldb), function(i) median(scaffoldb[nearest$nn.idx[i,c(1:maxk)]]))
    scaffold2
  }

  suppressMessages({
  suppressWarnings({
    seurat_obj <- Seurat::CreateSeuratObject(counts = as.matrix(ranked_df),  min.cells = 0, min.features = 0)
    seurat_obj <- Seurat::SetAssayData(object = seurat_obj, new.data = as.matrix(ranked_df))
  })
  Seurat::VariableFeatures(seurat_obj) <- top_genes
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = min(ncol(seurat_obj) - 1, max_pcs))

  threshold_var_explained <- 0.5
  stdev_explained <- Seurat::Stdev(object = seurat_obj, reduction = 'pca')
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


  return(predicted_df)

}

#' Binning
#'
#' This function performs a binning procedure by ranking cells within their predicted
#' potency category and arranging cells uniformly by rank per potency category
#' within equal length partitions of the unit interval, yielding binned smooth potency score.

#'
#' @param predicted_df dataframe of 4 columns: predicted potency scores and categories (CytoTRACE 2 model predicted values from predict function)
#' and smoothed potency score and potency categories (processed values from smoothData function)
#'
#' @return dataframe with updated predicted potency scores and categories
#' @export
#'
#' @examples
#' predicted_df --> output of smoothData
#' cytotrace <- binData(predicted_df)


binData <- function(predicted_df) {

  # message('cytotrace2: Started binning')
  
  labels <-
    c(
      'Differentiated',
      'Unipotent',
      'Oligopotent',
      'Multipotent',
      'Pluripotent',
      'Totipotent'
    )

  pred_potencies <- predicted_df$preKNN_CytoTRACE2_Potency
  unique_potency <- unique(pred_potencies)
  limits <- seq(0, 1, length.out = 7)
  score <- 'preKNN_CytoTRACE2_Score'
  for (potency_i in seq(1:length(labels))) {

    potency <- labels[potency_i]
    lower <- limits[potency_i]
    upper <- limits[potency_i + 1]

    if (potency %in% unique_potency) {
      indices <- which(predicted_df$preKNN_CytoTRACE2_Potency == potency)
      data_order <- predicted_df[indices, score]
      index <- indices[order(data_order)]
      n <- length(index)
      min_max_scaler <- function(x, feature_range = c((lower + 1e-8), (upper - 1e-8))) {
        if (length(x) == 1) {
          scaled_values <-  (feature_range[2] - feature_range[1])/2 + feature_range[1] #assign to the midpoint
        } else {
          scaled_values <- (x - min(x)) / (max(x) - min(x)) * (feature_range[2] - feature_range[1]) + feature_range[1]
        }
        return(scaled_values)
      }

      scaled <- min_max_scaler(1:n)
      predicted_df[index, score] <- scaled
    }
  }
  return(predicted_df)

}




#' Smoothing
#'
#' This function generates single-cell predictions of differentiation status from scRNA-seq data.
#' It takes in a matrix of gene expression values where columns are cells and rows are genes.\cr\cr
#' annotations matching the number of columns (cells) in the matrix. Furthermore, for the analysis of
#' large datasets (>3,000 cells), users can increase the speed performance by enabling a subsampling approach for calculation
#' and using multiple cores.
#'
#' @param mat matrix of ranked gene expression where rows are cells and columns are genes
#' @param predicted_df dataframe of 2 columns: predicted potency score 
#' and predicted potency category for each cell (CytoTRACE2 model predicted values from predict function, prior to postprocessing)
#' @param top_genes list of the top 1000 most variable genes (output of dispersion function)
#' @param parallelize_smoothing Logical, whether to run the function on subsamples in parallel on multiple threads (default is TRUE).
#' @param ncores integer indicating the number of cores to utilize
#' @param smooth_batch_size integer indicating the number of cells to subsample
#' @param seed Integer, specifying the seed for reproducibility in random processes (default is 14).

#' @return a numeric vector of the predicted potency score of single cells from 1.0 (least differentiated) to 0.0 (most differentiated)
#' @importFrom HiClimR fastCor
#' @importFrom parallel mclapply detectCores
#' @import doParallel
#' @export
#'
#' @examples
#' mat --> first element in the output of preprocessData()
#' predicted_df --> output of predictData()
#' SmoothScore <- smoothData(expression_data, predicted_df, ncores = 4, parallelize_smoothing = TRUE, smooth_batch_size = 1000)


smoothData <- function(mat, predicted_df, top_genes, parallelize_smoothing,
                       ncores, smooth_batch_size, seed){

  # message('cytotrace2: Started Postprocessing: Smoothing CytoTRACE 2 model predicted scores')

  if(nrow(mat) < 1000 && parallelize_smoothing == TRUE){
    parallelize_smoothing <- FALSE
    ncores <- 1
    message("The number of cells in your dataset is less than 1000. Fast mode has been disabled.")
  }

  if (parallelize_smoothing) {
    if(Sys.info()["sysname"] == "Windows") {
      ncores <- 1
      message("Windows OS can run only on 1 core")
    } else {
      if (is.null(ncores)) {
        ncores <- parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1
      }}
  } else {
    if (is.null(ncores)) {
      ncores <- 1
    }
  }

  # managing the smooth_batch_size
  if (smooth_batch_size > nrow(mat)) {
    message ( "The passed subsample size is greater than the number of cells in the subsample.\nNow setting subsample size to ", nrow(mat), "\nPlease consider reducing the smooth_batch_size to a number in range 1000 - 3000 for runtime and memory efficiency.")
    size <- nrow(mat)
  } else if  (nrow(mat) > 1000 && smooth_batch_size > 1000) {
    size <- smooth_batch_size
    message("Please consider reducing the smooth_batch_size to a number in range 1000 - 3000 for runtime and memory efficiency.")
  } else {
    size <- smooth_batch_size
  }

  # generating subsamples
  set.seed(seed)
  chunk_number <- round(nrow(mat)/size)
  subsamples <- split(seq_len(nrow(mat)), sample(factor(seq_len(nrow(mat)) %% chunk_number)))
 
  nn_pred_order <- predicted_df$preKNN_CytoTRACE2_Score
  names(nn_pred_order) <- rownames(predicted_df)
  sample_names <- lapply(subsamples, function(x) names(nn_pred_order)[x])

  sub_mat <- mat[, top_genes]


  get_score <- function(subsample) {
    mt <- sub_mat[subsample,]

    get_markov_matrix <- function(mat) {

      D <- HiClimR::fastCor(base::t(mat))  # Pairwise pearson-r corrs

      diag(D) <- 0
      D[is.na(D)] <- 0
      cutoff <- max(mean(D), 0)
      D[D < cutoff] <- 0

      A <- D / (rowSums(D) + 1e-5)
      return(A)
    }

    smoothing_by_diffusion <- function(score, markov_mat, maxiter = 1e4) {
      init_score <- score
      prev_score <- score

      for (i in 1:maxiter) {
        cur_score <- 0.9 * markov_mat %*% prev_score + 0.1 * init_score

        if (mean(abs(cur_score - prev_score)) / (mean(init_score) + 1e-6) < 1e-6) {
          break
        }

        prev_score <- cur_score
      }

      return(cur_score)
    }

    markov_matrix <- get_markov_matrix(mt)
    score <- smoothing_by_diffusion(nn_pred_order[subsample], markov_matrix)
    return(score)
  }

  if (parallelize_smoothing) {
    message('cytotrace2: Running with fast mode (subsamples are processed in parallel)')
    message(paste("This section will run on", chunk_number, "sub-sample(s) of approximately",
                  round(mean(unlist(lapply(subsamples, length)))), "cells each using", min(chunk_number, ncores),"/", parallel::detectCores(all.tests = FALSE, logical = TRUE), "core(s)."))
    cytotrace <- parallel::mclapply(subsamples, mc.cores = min(chunk_number, ncores), get_score)
  } else {
    message('cytotrace2: Running with slow mode (subsamples are processed sequentially)')

    cytotrace <- lapply(subsamples, get_score)
  }

  cytotrace <-  unlist(cytotrace)
  names(cytotrace) <- unlist(sample_names)


  #Final steps
  cytotrace <- cytotrace[rownames(predicted_df)]

  return(cytotrace)

}




