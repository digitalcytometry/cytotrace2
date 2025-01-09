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

#' Shortest Consensus Function
#' This function finds the smallest segment size in the predicted potency scores of the neighboring cells
#' where the average scores of two consecutive portions fall into the same potency category
#' and returns twice that segment length.
#' @param neighbor_scores A numeric vector of predicted potency scores of neighboring cells.
#' @return An integer representing  2*shortest segment size
#' @export

shortest_consensus <- function(neighbor_scores) {
  labels <- c('Differentiated', 'Unipotent', 'Oligopotent', 'Multipotent', 'Pluripotent', 'Totipotent')
  breaks <- seq(0, 1, length.out = 7) 

  len <- length(neighbor_scores)
  idx_use <- 2
  last_part <- FALSE
  max_i <- floor(len / 2)

  for (i in 2:max_i) {
    first_indices <- 1:i
    second_indices <- (i + 1):(2 * i)
        if (2 * i > len) {
      second_indices <- (i + 1):len
    }
    
    mean1 <- mean(neighbor_scores[first_indices], na.rm = TRUE)
    mean2 <- mean(neighbor_scores[second_indices], na.rm = TRUE)
    
    potency1 <- as.character(cut(mean1, 
                 breaks = breaks, 
                 labels = labels, 
                 include.lowest = TRUE, 
                 right = TRUE))
    potency2 <- as.character(cut(mean2, 
                 breaks = breaks, 
                 labels = labels, 
                 include.lowest = TRUE, 
                 right = TRUE))
    
    if (!is.na(potency1) && !is.na(potency2) &&
        potency1 == potency2 && !last_part) {
      idx_use <- i
      last_part <- TRUE
    }
  }
  
  return(2 * idx_use)
}

#' Scale function
#' This custom function scales the rows of a matrix to have a mean of 0 and a standard deviation of 1.
#' It handles cases where the standard deviation of a row is very small, and mean of the scaled row is not close to 0.
#' It uses population standard deviation instead of sample standard deviation.
#' @param x A matrix of numeric values where rows are cells and columns are genes.
#' @param constant_threshold A numeric value representing the threshold below which the standard deviation of a row is considered very small (default is 10 * .Machine$double.eps).
#' @param mean_tolerance A numeric value representing the tolerance for the mean of the scaled row to be considered close to 0 (default is 1e-10).
#' @return A matrix of scaled values where rows have a mean of 0 and a standard deviation of 1.
#' @export

scale_rows_custom <- function(x, constant_threshold = 10 * .Machine$double.eps, mean_tolerance = 1e-10) {
  row_means <- rowMeans(x, na.rm = TRUE) 
  
  ### changed to population standard deviation
  pop_sd <- function(x, na.rm = TRUE) {
    sqrt(mean((x - mean(x))^2, na.rm = na.rm))
  }
  row_sds <- apply(x, 1, pop_sd, na.rm = TRUE)
  
  # if sds are smaller than very small threshold, set to 1
  row_sds_handled <- ifelse(row_sds < constant_threshold, 1, row_sds)
  
  scaled_x <- sweep(x, 1, row_means, FUN = "-")
  scaled_x <- sweep(scaled_x, 1, row_sds_handled, FUN = "/")
  
  # mean_2 
  residual_means <- rowMeans(scaled_x, na.rm = TRUE)
  
  # rows where means of scaled rows are not close to 0
  rows_to_adjust <- which(abs(residual_means) > mean_tolerance)
  
  if (length(rows_to_adjust) > 0) {
    scaled_x[rows_to_adjust, ] <- sweep(scaled_x[rows_to_adjust, ], 1, residual_means[rows_to_adjust],FUN = "-")
  }
  
  return(scaled_x)
}
#' Smooth Data with kNN
#'
#' This applies adaptive nearest neighbor smoothing to datasets with >100 cells.
#' It identifies the nearest neighbors of each cell based on top principal components,
#' recalculates each cell’s potency score as a distance-weighted average of its neighbors 
#' including itself, and updates potency category accordingly.
#'
#' @param log2_df matrix of log2-transformed gene expression values where rows are cells and columns are genes
#' @param predicted_df dataframe of 4 columns containing CytoTRACE 2 model predicted and postprocessed potency scores and categories
#' @param seed integer specifying the seed for reproducibility in random processes (default is 14).
#'
#' @return A dataframe with final updated predictions
#' @importFrom Seurat CreateSeuratObject ScaleData RunPCA VariableFeatures Stdev
#' @importFrom RSpectra svds
#' @export
#'
#' @examples
#'  log2_data --> second element in the output of preprocessData()
#'  predicted_df --> output of binData()
#'  predicted_df <- smoothDatakNN(ranked_df, predicted_df, seed = 14, ncores = 4)

smoothDatakNN <- function(log2_df, predicted_df, seed, ncores){
  
  set.seed(seed)
  
  df_out <- copy(predicted_df)
  df_out$CytoTRACE2_Score <- 0.0
  df_out$CytoTRACE2_Potency <- ''
  
  message("Number of cores for KNN: ", ncores)
  
  # Perform truncated SVD
  cell_names_use <- rownames(log2_df)
  log2_data_scaled <- scale_rows_custom(log2_df)
  svd_res <- RSpectra::svds(log2_data_scaled, k = 30, opts = list(center = TRUE, scale = FALSE))
  df_pca <- svd_res$u %*% diag(svd_res$d)
  rownames(df_pca) <- cell_names_use
  
  # KNN smoothing over all cells
  num_cells <- length(cell_names_use) 
  chunk_size <- ceiling(num_cells/ncores)
  res_all <- parallel::mclapply(seq.int(1,ncores),mc.cores=ncores,function(chunk_id){
    idx_start <- (chunk_id-1)*chunk_size+1
    idx_end <- min(chunk_id*chunk_size,num_cells)
    cell_names_chunk <- cell_names_use[idx_start:idx_end]
    res_chunk <- sapply(cell_names_chunk, function(cell) {
      dist_col <- sqrt(rowSums(sweep(df_pca, 2, df_pca[cell,], "-")^2))
      dist_col <- dist_col / max(dist_col)
      top30_idx <- head(order(dist_col), 30)
      neighbor_dists <- dist_col[top30_idx]
      neighbor_scores <- df_out[names(neighbor_dists), 'preKNN_CytoTRACE2_Score']
      num_neighbors_keep <- shortest_consensus(neighbor_scores)
      
      if (num_neighbors_keep > 1) {
        keep_idx <- top30_idx[1:num_neighbors_keep]
        neighbor_dists <- dist_col[keep_idx]
        neighbor_scores <- df_out[names(neighbor_dists), 'preKNN_CytoTRACE2_Score']
        diff_vec <- (1 - neighbor_dists)^2
        proposed_new_score <- sum(neighbor_scores * diff_vec) / sum(diff_vec)
      } else {
        proposed_new_score <- df_out[cell, 'preKNN_CytoTRACE2_Score']
      }
      proposed_new_score
    })
    res_chunk
  })
  res <- Reduce(c,res_all)
    
  df_out[names(res), 'CytoTRACE2_Score'] <- res
  
  # Re-map potency categories based on new scores
  ranges <- seq(0, 1, length.out = 7)
  labels <- c('Differentiated', 'Unipotent', 'Oligopotent', 
              'Multipotent', 'Pluripotent', 'Totipotent')
  
  df_out$CytoTRACE2_Potency <- cut(df_out$CytoTRACE2_Score,
                                        breaks = ranges,
                                        labels = labels,
                                        include.lowest = TRUE,
                                        right = TRUE)
  
  return(df_out)
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


smoothData <- function(mat, predicted_df, top_genes, ncores, smooth_batch_size, parallelize_smoothing, seed){

  # message('cytotrace2: Started Postprocessing: Smoothing CytoTRACE 2 model predicted scores')

  if(nrow(mat) < 1000 && parallelize_smoothing == TRUE){
    parallelize_smoothing <- FALSE
    ncores <- 1
    message("The number of cells in your dataset is less than 1000. Fast mode has been disabled.")
  }

  # managing the smooth_batch_size
  if (smooth_batch_size > nrow(mat)) {
    message ( "The passed subsample size is greater than the number of cells in the subsample.\nNow setting subsample size to ", nrow(mat), "\nPlease consider reducing the smooth_batch_size to 1000 for runtime and memory efficiency.")
    size <- nrow(mat)
  } else if  (nrow(mat) > 1000 && smooth_batch_size > 1000) {
    size <- smooth_batch_size
    message("Please consider reducing the smooth_batch_size to 1000 for runtime and memory efficiency.")
  } else {
    size <- smooth_batch_size
  }

  # generating subsamples
  set.seed(seed)
  chunk_number <- ceiling(nrow(mat)/size)
  subsamples <- split(seq_len(nrow(mat)), sample(factor(seq_len(nrow(mat)) %% chunk_number)))
 
  nn_pred_order <- predicted_df$preKNN_CytoTRACE2_Score
  names(nn_pred_order) <- rownames(predicted_df)
  sample_names <- lapply(subsamples, function(x) names(nn_pred_order)[x])

  sub_mat <- mat[, top_genes]


  get_score <- function(subsample) {
    mt <- sub_mat[subsample,]

    get_markov_matrix <- function(mat) {

      D <- HiClimR::fastCor(base::t(mat), verbose = FALSE)  

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