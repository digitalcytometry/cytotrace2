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
#' @param full_model Logical, indicating whether to predict based on the full ensemble of 17 models
#' or a reduced ensemble of 5 most predictive models (default is FALSE).
#' @param parallelize_models Logical, indicating whether to run the prediction function on
#' models in parallel on multiple threads (default is TRUE).
#' @param ncores Integer, indicating the number of cores to utilize when parallelize_models
#' and/or parallelize_smoothing are TRUE (default is NULL, pipeline detects the number of
#' available cores and runs on that number; for Windows, it will be set to 1).
#' @param batch_size Integer or NULL, indicating the number of cells to subsample for the pipeline steps.
#'  No subsampling if NULL.(default is 10000, recommended value for input data size > 10K cells).
#' @param smooth_batch_size Integer or NULL, indicating the number of cells to subsample further
#' within the batch_size for the smoothing step of the pipeline. No subsampling if NULL.
#' (default is 1000, recommended value for input data size > 1K cells).
#' @param max_pcs Integer, indicating the maximum number of principal components to use in the smoothing by kNN step (default is 200).
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
#' Minji Kang, Jose Juan Almagro Armenteros, Gunsagar Gulati, Rachel Gleyzer, Erin Brown and Susanna Avagyan.
#'
#' @section 
#' [LICENSE]](https://github.com/digitalcytometry/cytotrace2/blob/master/LICENSE).
#' 
#' @export



cytotrace2 <- function(input,
                       species = "mouse",
                       is_seurat = FALSE,
                       slot_type = "counts",
                       full_model = FALSE,
                       batch_size = 10000,
                       smooth_batch_size = 1000,
                       parallelize_models = TRUE,
                       parallelize_smoothing = TRUE,
                       ncores = NULL,
                       max_pcs = 200,
                       seed = 14) {


  set.seed(seed)

  message('cytotrace2: Started loading data')
  
  # if input is Seurat object but the is_seurat is not set to TRUE
  if (class(input) == "Seurat" & is_seurat == FALSE) {
    stop("The input is a Seurat object. Please make sure to set is_seurat = TRUE.")
  }

  if (is.character(input) && file.exists(input)) { # check if the provided input is a filepath
    if (is_seurat == FALSE) { # check if the provided is a path to a seurat object
      data <- loadData(input)
    } else {
      seurat <- readRDS(input)
      data <- loadData_fromSeurat(seurat, slot_type)
    } } else {
      if (is_seurat == FALSE) { # check if the provided is a seurat object
        data <- copy(input)
      } else {
        seurat <- copy(input)
        data <- loadData_fromSeurat(input, slot_type)
      }
    }
  
  if (any(duplicated(colnames(data)))) {
    stop("Please make sure the cell/sample names are unique")
  }
  
  if (any(duplicated(rownames(data)))) {
    stop("Please make sure the gene names are unique")
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
        ncores <- parallel::detectCores(all.tests = FALSE, logical = TRUE) - 1
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
  chunk <- round(ncol(data) / batch_size)
  subsamples <- split(seq_len(ncol(data)), sample(factor(seq_len(ncol(data)) %% chunk)))
  sample_names <- lapply(subsamples, function(x) colnames(data)[x])

  message('cytotrace2: Running on ', chunk, " subsample(s) approximately of length " , batch_size)

  # reading preconstructed model weights
  if (full_model) {
    parameter_dict <-
      readRDS(system.file("extdata", "parameter_dict_17.rds",
                          package = "CytoTRACE2"))

  } else {
    parameter_dict <-
      readRDS(system.file("extdata", "parameter_dict_5_best.rds",
                          package = "CytoTRACE2"))
  }

  nc <- min(chunk, ncores)
  # composite function to preprocess and predict each subsample
  subsample_processing_f <-   function(subsample) {
    dt <- data[,subsample]
    message('cytotrace2: Started preprocessing.')
    # preprocessing
    ranked_data <- preprocessData(dt, species)
    gene_names <- colnames(ranked_data)
    cell_names <- rownames(ranked_data)
    
    message('cytotrace2: Started prediction.')
    # predicting
    predicted_df <-
      predictData(parameter_dict, ranked_data, parallelize_models, ncores = nc)
    rownames(predicted_df) <- cell_names

    # calculating top 1000 most variable genes
    num_genes <- ncol(ranked_data)
    dispersion_index <- sapply(1:num_genes, function(i) disp_fn(ranked_data[, i]))
    top_genes <- gene_names[order(dispersion_index, decreasing = TRUE)[1:min(1000, num_genes)]]
    # top_genes <- colnames(ranked_data)[top_genes_idx]
    
    message('cytotrace2: Started postprocessing.')
    # smoothing
    smoothScore <-
      smoothData(
        ranked_data,
        predicted_df,
        top_genes,
        ncores = ncores,
        smooth_batch_size = smooth_batch_size,
        parallelize_smoothing = parallelize_smoothing,
        seed = seed
      )

    predicted_df$preKNN_CytoTRACE2_Score <- smoothScore

    # running postprocessing
    if (nrow(ranked_data) <= 10) {
      message('cytotrace2: Number of cells in data is less than 10. Skipping postprocessing.')
      predicted_df$CytoTRACE2_Potency <- predicted_df$preKNN_CytoTRACE2_Potency
      predicted_df$CytoTRACE2_Score <- predicted_df$preKNN_CytoTRACE2_Score
      
    } else {

      ## binning method to fit smoothed values back to original predictions
      predicted_df <- binData(predicted_df)
      predicted_df <- predicted_df[cell_names,]


      if (nrow(ranked_data) <= 100)  {
        message('cytotrace2: Number of cells in data is less than 100. Skipping kNN smoothing.')
        predicted_df$CytoTRACE2_Potency <- predicted_df$preKNN_CytoTRACE2_Potency
        predicted_df$CytoTRACE2_Score <- predicted_df$preKNN_CytoTRACE2_Score
      }
      if (sd(ranked_data) == 0) {
        message('cytotrace2: Zero variance of ranked matrix. Skipping kNN smoothing.')
        predicted_df$CytoTRACE2_Potency <- predicted_df$preKNN_CytoTRACE2_Potency
        predicted_df$CytoTRACE2_Score <- predicted_df$preKNN_CytoTRACE2_Score
      }
      if (nrow(ranked_data) > 100 && sd(ranked_data) != 0) {
        # message('cytotrace2: Started postprocessing: Smoothing by kNN')
        ranked_df <- data.frame(base::t(ranked_data))
        rownames(ranked_df) <- gene_names
        colnames(ranked_df) <- cell_names
        predicted_df <- smoothDatakNN(ranked_df,
                                      predicted_df,
                                      top_genes,
                                      max_pcs,
                                      seed)
        predicted_df <- predicted_df[cell_names,]
      }
    }


    predicted_df

  }

  message('cytotrace2: Started running on subsample(s). This will take a few minutes.')

  results <- lapply(subsamples,subsample_processing_f)

  predicted_df <- do.call(rbind, results)
  rownames(predicted_df) <- unlist(sample_names)
  predicted_df <- predicted_df[init_order,]

  # add relative score
  ranked_scores <- rank(predicted_df$CytoTRACE2_Score)
  predicted_df$CytoTRACE2_Relative <- (ranked_scores - min(ranked_scores)) / (max(ranked_scores) - min(ranked_scores))
  
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

