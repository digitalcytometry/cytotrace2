#' Binary Module
#'
#' This function generates predictions for each layer using the input parameter dictionary for a given layer
#'
#' @param layer_dict The layer dictionary containing weights and parameters.
#' @param ranked_data The input ranked data.
#' @param eps The epsilon value for numerical stability in calculations (default: 1e-05).
#' @return A list containing the calculated matrix R_all before batch normalization and the vector pred of predictions
#' @importFrom Matrix Matrix tcrossprod rowSums
#' @export
#'
#' @examples
#' output <- BinaryModule(layer_dict, ranked_data)

BinaryModule <- function(layer_dict, ranked_data, eps = 1e-05) {

  alpha <- layer_dict$weight
  n <- layer_dict$n
  maxrank <- layer_dict$maxrank_norm

  x <- pmin(ranked_data, maxrank)
  x_counts <- x < maxrank

  x <- Matrix::Matrix(x, sparse = TRUE)
  x_counts <- Matrix::Matrix(x_counts, sparse = TRUE)
  alpha <- Matrix::Matrix(alpha, sparse = TRUE)

  R <- Matrix::tcrossprod(x, Matrix::t(alpha))
  R_counts <- Matrix::tcrossprod(x_counts, Matrix::t(alpha)) / (Matrix::tcrossprod(matrix(Matrix::rowSums(x_counts), ncol = 1), base::t(n/dim(x_counts)[2])))

  R_out <- 1 - sweep(sweep(R,2,c((n*(n+1))/2)), 2, (c(n) * maxrank), "/")
  R_all <- cbind(R_out, R_counts)

  R_out_norm <- sweep(sweep(R_all,2, layer_dict$running_mean, '-'), 2, sqrt(layer_dict$running_var + eps), "/" )
  pred <- Matrix::tcrossprod(layer_dict$out.weight, R_out_norm) + layer_dict$out.bias

  return(pred@x)

}


#' Softmax
#'
#' This function applies the softmax function to the input vector.
#'
#' @param x The input vector.
#' @return The vector after applying softmax.
#' @export
#'
#' @examples
#' softmax_values <- softmax(x)

softmax <- function(x) {
  return(exp(x) / sum(exp(x)))
}

#' Binary Encoder
#'
#' This function generates predicted probabilities and predicted potency score of each sample based on a given model
#' by iterating over all layers of that model and aggregating results
#'
#' @param model_dict The model dictionary containing layer_dictionaries of weights and parameters.
#' @param ranked_data The input ranked data.
#' @return A list containing the sample_size X layer_number sized prob_pred matrix of predicted probabilities and
#' sample_size sized vector of predicted potency score containing values in the range 0-1.
#' @export
#'
#' @examples
#' encoded_output <- BinaryEncoder(model_dict, ranked_data)

BinaryEncoder <- function(model_dict, ranked_data) {

    pred_list <- list()

  pred_list <- sapply(1:6, function(layer) {
    layer_dict <- model_dict[[layer]]
    BM_output <- BinaryModule(layer_dict, ranked_data)
    BM_output
  }, simplify = "list")

  prob_pred <- base::t(apply(pred_list, 1, softmax))

  order_vector <- matrix(0:5 / 5, 6, nrow = 6, ncol = 1)
  prob_order <- as.matrix(prob_pred) %*% order_vector

  return(list('prob_pred' = prob_pred, 'prob_order' = prob_order))

}


#' Prediction
#'
#' This function generates predicted probabilities and predicted potency score of each sample based on a given
#' ensemble of models by iterating over all models of the set of models and aggregating results.
#' It then performs post processing of predicted probabilities and predicted potency scores
#' to generate 1 categorical label and 1 continuous score value for each sample
#' by aggregating results over all models and layers.
#'
#' @param parameter_dict The ensemble parameter dictionary containing model_dictionaries of weights and parameters.
#' @param ranked_data The input ranked data.
#' @param parallelize_models boolean value, if TRUE models of the ensemble will run on multiple threads
#' @param ncores integer indicating the number of cores to utilize
#'
#' @return The dataframe containing predicted potency category and score for each sample (with sample names as rownames) prior to postprocessing
#' @importFrom parallel mclapply detectCores
#' @importFrom dplyr bind_cols
#' @import doParallel
#' @export
#'
#' @examples
#' parameter_dict <- readRDS(system.file("extdata", "parameter_dict_17.rds", package = "CytoTRACE2"))
#' predicted_df <- predictData(parameter_dict, ranked_data, parallelize_models = TRUE, ncores = 4)

predictData <- function(parameter_dict, ranked_data, parallelize_models, ncores) {

  # message('cytotrace2: Started prediction')

  if (parallelize_models) {
    message(paste("This section will run using ", ncores,"/", parallel::detectCores(all.tests = FALSE, logical = TRUE), "core(s)."))
    results <- parallel::mclapply(parameter_dict, mc.cores = ncores, function(model_dict) {
      BE_output <- BinaryEncoder(model_dict, ranked_data)
      cbind(BE_output$prob_pred,
            BE_output$prob_order)
    })

  } else {

    results <- lapply(parameter_dict, function(model_dict) {
      BE_output <- BinaryEncoder(model_dict, ranked_data)
      cbind(BE_output$prob_pred,
            BE_output$prob_order)
    })
  }

  results_binded <- as.matrix(dplyr::bind_cols(results))
  is_order_col <- seq_len(ncol(results_binded)) %% 7 == 0

  all_order_out <- results_binded[, is_order_col, drop = FALSE]
  all_preds_out <- results_binded[, !is_order_col, drop = FALSE]

  predicted_score <- rowMeans(all_order_out)


  all_preds_out_avg <- list()
  for (i in 0:5) {
    col_indices <- seq(from = 1+i, to = ncol(all_preds_out), by = 6)
    row_means <- rowMeans(as.matrix(all_preds_out[, col_indices]))
    all_preds_out_avg <- cbind(all_preds_out_avg, row_means)
  }

  colnames(all_preds_out_avg) <- seq(1:6)

  predicted_category <- max.col(all_preds_out_avg)

  labels <- c('Differentiated', 'Unipotent', 'Oligopotent', 'Multipotent', 'Pluripotent', 'Totipotent')
  labels_dict <- setNames(labels, 1:6)

  predicted_df <- data.frame('preKNN_CytoTRACE2_Score' = predicted_score, 'preKNN_CytoTRACE2_Potency' = predicted_category)
  predicted_df$preKNN_CytoTRACE2_Potency <- labels_dict[predicted_df$preKNN_CytoTRACE2_Potency]
  row.names(predicted_df) <- rownames(ranked_data)

  return(predicted_df)
}

