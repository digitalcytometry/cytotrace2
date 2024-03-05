#' Load Data from File
#'
#' This function loads data from a specified file path.
#'
#' @param filepath The path to the file containing the data
#'    (count matrix, samples as columns, genes as rows)
#' @return The loaded data as a data frame.
#' @importFrom data.table fread
#' @export
#'
#' @examples
#' data <- loadData("/path/to/file.tsv")
#'

loadData <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("File not found.")
  }

  # load data
  data <- data.table::fread(filepath, header = TRUE,  sep = "\t", check.names = FALSE, data.table = FALSE)

  # check for uniqueness of cell and gene names
  if (any(duplicated(data[[1]]))) {
    stop("Please make sure the gene names are unique")
  }
  
  row.names(data) <- data[[1]]
  row.names(data) <- gsub(".", "-", row.names(data), fixed = T)
  data[[1]] <- NULL
  
  if (any(duplicated(colnames(data)))) {
    stop("Please make sure the cell/sample names are unique")
  }
 

  return(data)
}


#' Load Data from Seurat object
#'
#' This function loads data from a Seurat object
#'
#' @param object The Seurat object
#' @param slot_type The type of slot to access from "RNA" assay. Can be 'counts' or 'data' (default is 'counts')
#' @return The RNA assay data as a data frame.
#' @importFrom Seurat GetAssayData
#' @export
#'
#' @examples
#' devtools::install_github('satijalab/seurat-data')
#' library(SeuratData)
#' InstallData("pbmc3k")
#' data("pbmc3k")
#' data <- loadData_fromSeurat(pbmc3k)

loadData_fromSeurat <- function(object, slot_type) {

  # load data
  data <- as.data.frame(Seurat::GetAssayData(object = object, assay="RNA", slot=slot_type))

  return(data)
}



#' Preprocess Data
#'
#' This function maps the gene names (rownames) of input data to mouse gene names
#' based on a provided mapping dictionary, filters and orders the genes based on provided feature set,
#' and ranks the genes for each sample.
#' Caution: Please make sure that the input dataset contains gene symbols from only one species (either human or mouse),
#' following the gene symbol convention from HGNC for human and MGI for mouse.
#'
#' @param data The input data to be preprocessed, output of loadData function
#'   (count matrix, samples as columns, genes as rows)
#' @param species String, indicates the species name from which the gene names for the input data come from
#' (2 options: human or mouse)
#' @return Mapped and ranked expression matrix (samples as rows, genes as columns)
#' @importFrom data.table fread setDT setDF setnafill setkey
#' @importFrom plyr mapvalues
#' @importFrom Matrix t
#' @importFrom Rfast colRanks
#' @export
#'
#' @examples
#' ranked_data <- preprocessData(dt)
#'

preprocessData <- function (data, species) {
  
  # check for uniqueness of cell and gene names
  if (any(duplicated(data$V1))) {
    stop("Please make sure the gene names are unique")
  }
  
  if (any(duplicated(colnames(data)))) {
    stop("Please make sure the cell/sample names are unique")
  }

  if (!is.data.frame(data) & !is.data.table(data)) {
   message("The function expects an input of type 'data.frame' or 'data.table'.\nAttempting to convert the provided input to the required format.")
   data <- as.data.frame(data)
  }
  
  
  # Load the features_model.csv file
  features <- read.csv(system.file("extdata", "features_model_training_17.csv",
                                   package = "CytoTRACE2"),  row.names = 1,
                                   check.names = FALSE)[[1]]

  # mapping
  gene_names <- rownames(data)
  expression <- copy(data)
  data.table::setDT(expression)



  if (species == "human") {
    # Load the presaved mouse-human orthology table
    mt_dict <- data.table::fread(system.file("extdata", "mt_dict_human_to_mouse.csv",
                                             package = "CytoTRACE2"), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    # message("cytotrace2: Mapping input gene names to mouse orthologs")
    mapping <- unlist(plyr::mapvalues(gene_names, colnames(mt_dict), mt_dict[1,], warn_missing = FALSE))

    # map those unmapped genes that correspond to aliases of a human gene that has orthology but isn't in the dataset
    rescue_dict <- read.csv(system.file("extdata", "mt_alias_dict.csv",
                                        package = "CytoTRACE2"), header = TRUE, check.names = FALSE)
    map_df <- data.frame("original_gene" = gene_names, 'mapped_gene' = mapping, "mapped" = ifelse(mapping %in% mt_dict[1,], "1", "0"))
    mapped_genes <- map_df[map_df$mapped == 1,]$original_gene
    unmapped_genes <- map_df[map_df$mapped == 0,]$original_gene
    is_alias <- intersect(unmapped_genes, rescue_dict$alias)
    original_to_alias <- rescue_dict[rescue_dict$alias %in% is_alias,]$hsgene
    mapped_original_to_alias <- intersect(original_to_alias, mapped_genes)
    alias_to_unmapped_original <- rescue_dict[rescue_dict$hsgene %in% setdiff(original_to_alias, mapped_original_to_alias),]$alias
    rownames(rescue_dict) <- rescue_dict$alias
    unmapped_original_with_alias <- map_df[map_df$original_gene %in% alias_to_unmapped_original,]$mapped_gene
    map_df[map_df$original_gene %in% alias_to_unmapped_original,]$mapped_gene <- rescue_dict[unmapped_original_with_alias,]$mmgene
    mapping <- map_df$mapped_gene

  } else {
    mapping <- gene_names
  }


  expression[, mapped_genes:= mapping]
  data.table::setkey(expression, mapped_genes)
  
  
  # message("cytotrace2: Reducing genes to preselected feature set")
  message(length(base::intersect(mapping, features)), " input genes mapped to model genes.")
  
  # if mapped genes are too few, warn for possible incompatible dataset or wrong species input
  if (length(base::intersect(mapping, features)) < 9000) {
    warning("The number of input genes mapped to the model is too low. Please verify the input species is correct.\nIn case of a correct species input, be advised that model performance might be compromised due to gene space differences.")
  }
  # reducing to pre-selected features
  expression_final <- expression[.(features) ]
  data.table::setnafill(expression_final, fill = 0, cols=colnames(expression_final)[-ncol(expression_final)])

  # message("cytotrace2: Ranking gene expression values within each cell/sample.")

  data.table::setDF(expression_final)
  gene_names <- expression_final$mapped_genes
  expression_final$mapped_genes <- NULL
  cell_names <- colnames(expression_final)


  ranked_data <- base::t(Rfast::colRanks(as.matrix(expression_final), descending = TRUE, method = "average"))
  colnames(ranked_data) <- gene_names
  rownames(ranked_data) <- cell_names

  # return(list(expression_final, ranked_data))

  return(ranked_data)


}
