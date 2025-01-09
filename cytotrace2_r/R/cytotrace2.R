#' @title cytotrace2: Predicting Cellular Potency Score and Category
#'
#' @description This function generates single-cell predictions of developmental potential from scRNA-seq data.
#' It takes in a matrix of gene expression values where columns are cells and rows are genes \cr\cr
#'
#' @param input The input single-cell RNA-seq data. It can be either an expression matrix
#' (rows as genes, columns as cells) or the filepath to .tsv file containing the data table,
#' or a Seurat object or the filepath to an .rds file containing a Seurat object.
#' @param species String, indicates the species name from which the gene names for the
#' input data come from (two options: "human" or "mouse", default is "mouse").
#' @param is_seurat Logical, indicating whether the input is a Seurat object/filepath to a
#' Seurat object (default is FALSE).
#' @param slot_type Character, indicating the type of slot to access from "RNA" assay if provided is a Seurat object &
#' is_seurat = TRUE. Can be 'counts' or 'data' (default is 'counts')
#' @param parallelize_smoothing Logical, indicating whether to run the smoothing function
#' on subsamples in parallel on multiple threads (default is TRUE).
#' @param parallelize_models Logical, indicating whether to run the prediction function on
#' models in parallel on multiple threads (default is TRUE).
#' @param ncores Integer, indicating the number of cores to utilize when parallelize_models
#' and/or parallelize_smoothing are TRUE (default is NULL, pipeline detects the number of
#' available cores and runs on half of them; for Windows, it will be set to 1).
#' @param batch_size Integer or NULL, indicating the number of cells to processat once, including subsampling 
#' for KNN smoothing. No subsampling if NULL.(default is 10000, recommended value for input data size > 10K cells).
#' @param smooth_batch_size Integer or NULL, indicating the number of cells to subsample further
#' within the batch_size for the smoothing by diffusion step of the pipeline. No subsampling if NULL.
#' (default is 1000, recommended value for input data size > 1K cells).
#' @param seed Integer, specifying the seed for reproducibility in random processes (default is 14).
#'
#' @return If `is_seurat` is FALSE, the function returns a dataframe containing predicted
#' raw and postprocessed cellular potency category and potency scores
#' If `is_seurat` is TRUE, adds metadata elements for all columns of prediction dataframe,
#' and returns the modified Seurat object.
#'
#' @details
#' The cytotrace2 function performs the following steps:
#'
#' 1. Loads and preprocesses the input data.
#' 2. Predicts single-cell developmental potential (with discrete potency categories and continuous (range 0-1) score).
#' 3. Smoothes the predicted values using diffusion and rescales using a combination of binning and smoothing by kNN approaches.
#'
#'
#' @examples
#' # example run with a scRNA-seq dataset (10x) encompassing 2280 cells from murine pancreatic epithelium (Bastidas-Ponce et al., 2019)
#' # download the .rds file (this will download the file to your working directory)
#' download.file("https://drive.google.com/uc?export=download&id=1ivi9TBlmzVTDGzNWQrXXeyL68Wug989K", "Pancreas_10x_downsampled.rds")
#' # load rds
#' data <- readRDS("Pancreas_10x_downsampled.rds")
#' # extract expression data
#' expression_data <- data$expression_data
#' # running CytoTRACE 2 main function - cytotrace2 - with default parameters
#' cytotrace2_result <- cytotrace2(expression_data)
#' 
#' @import magrittr
#' @import dplyr
#' @import stringr
#' @import data.table
#' @importFrom parallel detectCores
#' @importFrom Seurat AddMetaData
#' @importFrom Matrix t
#'
#' @seealso
#' See https://github.com/digitalcytometry/cytotrace2 for a detailed package description and instructions 
#'
#' @references
#' Link to the publication will be added upon being published
#'
#' @author
#' Minji Kang, Erin Brown, Jose Juan Almagro Armenteros, Gunsagar Gulati, Rachel Gleyzer, and Susanna Avagyan.
#'
#' @section 
#' [LICENSE]](https://github.com/digitalcytometry/cytotrace2/blob/master/LICENSE).
#' 
#' @export



cytotrace2 <- function(input,
                       species = "mouse",
                       is_seurat = FALSE,
                       slot_type = "counts",
                       batch_size = 10000,
                       smooth_batch_size = 1000,
                       parallelize_models = TRUE,
                       parallelize_smoothing = TRUE,
                       ncores = NULL,
                       seed = 14) {


  set.seed(seed)

  message('cytotrace2: Started loading data')

  # if input is Seurat object but the is_seurat is not set to TRUE
  if (class(input)[1] == "Seurat" & is_seurat == FALSE) {
    stop("The input is a Seurat object. Please make sure to set is_seurat = TRUE.")
  }

  if (is.character(input) && file.exists(input)) { # check if the provided input is a filepath
    if (is_seurat == FALSE) { # check if the provided is a path to a seurat object
      data <- loadData(input)
    } else {
      seurat <- readRDS(input)
      data <- loadData_fromSeurat(seurat, slot_type)
    } 
  } else {
    if (is_seurat == FALSE) { # check if the provided is a seurat object
      data <- copy(input)
    } else {
      seurat <- copy(input)
      data <- loadData_fromSeurat(input, slot_type)
    }
  }
  
  # check if the input corresponds to expected criteria
  # check for uniqueness of cell and gene names
  if (any(duplicated(rownames(data)))) {
    stop("Please make sure the gene names are unique")
  }
  
  if (any(duplicated(colnames(data)))) {
    stop("Please make sure the cell/sample names are unique")
  }

  # check for the format of the input
  if (!is.data.frame(data)) {
   message("The function expects an input of type 'data.frame' with gene names as row names and cell IDs as column names.\nAttempting to convert the provided input to the required format.")
   data <- as.data.frame(data)
  }

  # check if rownames are gene names (numbers as characters are not allowed)
  if (!all(grepl("[a-zA-Z]", rownames(data)))) {
    warning("The rownames of the input data are numeric. Please make sure the rownames are gene names.")
  }
  
  message("Dataset contains ", dim(data)[1], " genes and ", dim(data)[2], " cells.")
  
  # check if species are passed correctly
  is_human <- sum(sapply(rownames(data), function(x) all(toupper(x) == x))) / nrow(data) > 0.9
  is_mouse <- sum(sapply(rownames(data), function(x) all(toupper(x) != x))) / nrow(data) > 0.9
  
  if (is_human & species == 'mouse') {
    warning("Species is most likely human. Please revise the 'species' input to the function.")
  }
  
  if (is_mouse & species == 'human') {
    warning("Species is most likely mouse. Please revise the 'species' input to the function.")
  }
  
  if (max(data) < 20) {
      warning("It looks like your data may already be log-transformed. Please provide an untransformed expression matrix of raw counts or CPM/TPM normalized counts.")
  }

  # if user doesn't want to split into batches
  if (is.null(batch_size)) {
    batch_size <- ncol(data)
  }
  if (is.null(smooth_batch_size)) {
    smooth_batch_size <- ncol(data)
    parallelize_smoothing == FALSE
  }

  # if parallelization is enabled for data with <= 1000 cells
  if(ncol(data) <= 1000 && parallelize_smoothing == TRUE){
    parallelize_smoothing <- FALSE
    message("The number of cells in your dataset is less than 1000. Fast mode has been disabled.")
  }

  # managing the number of cores to be utilized
  if (parallelize_smoothing || parallelize_models) {
    if (Sys.info()["sysname"] == "Windows") {
      ncores = 1
      message("Windows OS can run only on 1 core")
    } else {
      if (is.null(ncores)) {
        ncores <- max(1, parallel::detectCores(all.tests = FALSE, logical = TRUE) %/% 2)
      }
    }
  } else {
    if (is.null(ncores)) {
      ncores <- 1
    }
  }



  # managing the subsamplesize for first block
  if (batch_size > ncol(data)) {
    message ("The passed subsample size is greater than the number of cells in dataset.\nNow setting subsample size to ", ncol(data))
    batch_size <- ncol(data)
  } else if  (ncol(data) > 10000 && batch_size > 10000) {
    message("Please consider reducing the batch_size to 10000 for runtime and memory efficiency.")
  }
  
  
  # fix the order of the input
  init_order <- colnames(data)

  # generating subsamples
  chunk <- ceiling(ncol(data) / batch_size)
  if (chunk > 1) {
    subsamples <- split(seq_len(ncol(data)), sample(factor(seq_len(ncol(data)) %% chunk)))
  } else {
    # No shuffling: return all columns in one group, preserving order
    subsamples <- list(seq_len(ncol(data)))
  }
  sample_names <- lapply(subsamples, function(x) colnames(data)[x])

  message('cytotrace2: Running on ', chunk, " subsample(s) approximately of length " , batch_size)

  parameter_dict <- readRDS(system.file("extdata", "parameter_dict_19.rds", package = "CytoTRACE2"))
  # B <- readRDS(system.file("extdata", "B_background.rds", package = "CytoTRACE2"))

  nc <- ncores

  # composite function to preprocess and predict each subsample
  subsample_processing_f <-  function(subsample) {
    dt <- data[,subsample]
    message('cytotrace2: Started preprocessing.')
    # preprocessing
    input_data <- preprocessData(dt, species)
    ranked_data <- input_data[[1]]
    log2_data <- input_data[[2]] 
    count_cells_few_genes <- input_data[[3]]

    gene_names <- colnames(ranked_data)
    cell_names <- rownames(ranked_data)
    message('cytotrace2: Started prediction.')
    # predicting
    predicted_df <-
      predictData(parameter_dict, ranked_data, log2_data, parallelize_models, ncores = nc)

    predicted_df <- predicted_df[cell_names,]
    # calculating top 1000 most variable genes
    num_genes <- ncol(log2_data)
    dispersion_index <- sapply(1:num_genes, function(i) disp_fn(log2_data[, i]))
    top_genes <- gene_names[order(dispersion_index, decreasing = TRUE)[1:min(1000, num_genes)]]
    # top_genes <- colnames(ranked_data)[top_genes_idx]
    
    message('cytotrace2: Started postprocessing.')
    # smoothing
    smoothScore <-
      smoothData(
        log2_data,
        predicted_df,
        top_genes,
        ncores = ncores,
        smooth_batch_size = smooth_batch_size,
        parallelize_smoothing = parallelize_smoothing,
        seed = seed
      )

    smoothScore <- smoothScore[cell_names]
    predicted_df$preKNN_CytoTRACE2_Score <- smoothScore

    # running postprocessing
    if (nrow(log2_data) <= 10) {
      message('cytotrace2: Number of cells in data is less than 10. Skipping postprocessing.')
      predicted_df$CytoTRACE2_Potency <- predicted_df$preKNN_CytoTRACE2_Potency
      predicted_df$CytoTRACE2_Score <- predicted_df$preKNN_CytoTRACE2_Score
      
    } else {

      ## binning method to fit smoothed values back to original predictions
      predicted_df <- binData(predicted_df)
      predicted_df <- predicted_df[cell_names,]

      if (nrow(log2_data) <= 100)  {
        message('cytotrace2: Number of cells in data is less than 100. Skipping kNN smoothing.')
        predicted_df$CytoTRACE2_Potency <- predicted_df$preKNN_CytoTRACE2_Potency
        predicted_df$CytoTRACE2_Score <- predicted_df$preKNN_CytoTRACE2_Score
      }
      if (sd(log2_data) == 0) {
        message('cytotrace2: Zero variance of ranked matrix. Skipping kNN smoothing.')
        predicted_df$CytoTRACE2_Potency <- predicted_df$preKNN_CytoTRACE2_Potency
        predicted_df$CytoTRACE2_Score <- predicted_df$preKNN_CytoTRACE2_Score
      }
      if (nrow(log2_data) > 100 && sd(log2_data) != 0) {
        # message('cytotrace2: Started postprocessing: Smoothing by kNN')

        predicted_df <- smoothDatakNN(log2_data,
                                      predicted_df,
                                      seed,
                                      ncores)
        predicted_df <- predicted_df[cell_names,]
      }
    }

    list(
      predicted_df = predicted_df,
      count_cells_few_genes = count_cells_few_genes
    )

  }

  message('cytotrace2: Started running on subsample(s). This will take a few minutes.')

  results <- lapply(subsamples,subsample_processing_f)
  
  all_counts <- sapply(results, `[[`, "count_cells_few_genes")
  all_count_cells_few_genes <- sum(all_counts)
  frac_cells_few_genes <- all_count_cells_few_genes / ncol(data)
  if (frac_cells_few_genes >= 0.2) {
    warning(sprintf(
      "WARNING: %.2f%% of input cells express fewer than %d genes. 
    For best results, a minimum gene count of 500-1000 is recommended. 
    Please see FAQ for guidelines at https://github.com/digitalcytometry/cytotrace2#frequently-asked-questions",
      frac_cells_few_genes * 100, 500
    ))
  }
  
  predicted_dfs <- lapply(results, `[[`, "predicted_df")
  predicted_df <- do.call(rbind, predicted_dfs)
  rownames(predicted_df) <- unlist(sample_names)
  predicted_df <- predicted_df[init_order,]

  # add relative score
  # ranked_scores <- rank(predicted_df$CytoTRACE2_Score)
  predicted_df$CytoTRACE2_Relative <- (predicted_df$CytoTRACE2_Score - min(predicted_df$CytoTRACE2_Score)) / (max(predicted_df$CytoTRACE2_Score) - min(predicted_df$CytoTRACE2_Score))
  predicted_df <- predicted_df[c("CytoTRACE2_Score", "CytoTRACE2_Potency" , "CytoTRACE2_Relative", "preKNN_CytoTRACE2_Score",  "preKNN_CytoTRACE2_Potency")]
  
  
  if (!is_seurat) {

    message('cytotrace2: Finished')

    return(predicted_df)

  } else {
    
    predicted_df <- predicted_df[colnames(seurat), ]
    seurat <- AddMetaData(object = seurat,
                          metadata = predicted_df)
    

    message('cytotrace2: Finished')

    return(seurat)
  }

}