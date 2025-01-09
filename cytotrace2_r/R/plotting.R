#' Plot Data
#'
#' This function generates plots (UMAPs, Boxplots) based on the input expression,
#' prediction, and annotation data by functions of the Seurat library
#'
#' @param cytotrace2_result output of cytotrace2 function (a dataframe with genes as rows and cells as columns,
#' or a Seurat object containing predictions of raw and final potency scores and categories)
#' @param expression_data Dataframe, expression data to be used for plotting
#' If cytotrace2 is a Seurat object and is_seurat = TRUE, can be left NULL (default is NULL).
#' @param annotation Dataframe, contains annotation data ' (optional, default is NULL, used only in phenotype based plots)
#'                  Has rownames which are single cell IDs as in expression data, and at least 1 columns for the phenotype/grouping variable of interest
#' @param is_seurat Logical, indicating whether the input is a Seurat object (default is FALSE). 
#' If TRUE and Seurat object contains PCA and UMAP embeddings, those will be automatically used for visualization, otherwise, PCA and UMAP will be calculated.
#' @param slot_type Character, indicating the type of slot to access from "RNA" assay if provided is a Seurat object &
#' is_seurat = TRUE. Can be 'counts' or 'data' (default is 'counts')
#' @param pc_dims The number of principal components to use for UMAP visualization (default: 30).
#' @param seed integer specifying the seed for reproducibility in random processes (default is 14).

#' @return A named list of 3 or 5 plots (3 UMAPs colored by predicted potency group, potency score and relative order,
#' and, if annotation file provided, 1 UMAP colored by phenotype and 1 boxplot of phenotype ~ potency score).
#' @import Seurat
#' @import ggplot2
#' @import dplyr
#' @export
#'
#' @examples
#'  # download the .rds file (this will download the file to your working directory)
#'  download.file("https://drive.google.com/uc?export=download&id=1ivi9TBlmzVTDGzNWQrXXeyL68Wug989K", "Pancreas_10x_downsampled.rds")
#'  # load rds
#'  data <- readRDS("Pancreas_10x_downsampled.rds")
#'  # extract expression data
#'  expression_data <- data$expression_data
#'  # running CytoTRACE 2 main function - cytotrace2 - with default parameters
#'  cytotrace2_result <- cytotrace2(expression_data)
#'  # extract annotation data
#'  annotation <- data$annotation
#'  # generate prediction and phenotype association plots with plotData function
#'  plots <- plotData(cytotrace2_result, expression_data, annotation)
#'  plots  
#'
plotData <- function(cytotrace2_result,
                     annotation = NULL,
                     expression_data = NULL,
                     is_seurat = FALSE,
                     slot_type = "counts",
                     pc_dims = 30,
                     seed = 14) {
  
  set.seed(seed)
  
  # Check if input is Seurat object but is_seurat is not set to TRUE
  if (class(cytotrace2_result) == "Seurat" & !is_seurat) {
    stop("The input is a Seurat object. Please make sure to set is_seurat = TRUE.")
  }
  
  # Check if expression data is not passed with non-Seurat input
  if (!is_seurat & is.null(expression_data)) {
    message("Please provide the expression data.\nIf using a Seurat object, set is_seurat = TRUE.")
  }
  
  message("Preparing input for visualization.")
  suppressMessages({
    suppressWarnings({
      
      if (is_seurat) {
        prediction <- data.frame(
          CytoTRACE2_Score = cytotrace2_result@meta.data$CytoTRACE2_Score,
          CytoTRACE2_Potency = cytotrace2_result@meta.data$CytoTRACE2_Potency,
          CytoTRACE2_Relative = cytotrace2_result@meta.data$CytoTRACE2_Relative,
          preKNN_CytoTRACE2_Potency = cytotrace2_result@meta.data$preKNN_CytoTRACE2_Potency,
          preKNN_CytoTRACE2_Score = cytotrace2_result@meta.data$preKNN_CytoTRACE2_Score
        )
        rownames(prediction) <- colnames(cytotrace2_result)
        
        if ("umap" %in% names(cytotrace2_result@reductions) & "pca" %in% names(cytotrace2_result@reductions)) {
          message("UMAP and PCA slots are found in the input Seurat object. Using them for visualization.")
          seurat <- copy(cytotrace2_result)
        } else {
          expression_data <- as.data.frame(GetAssayData(cytotrace2_result, assay = "RNA", slot = slot_type))
          
          if (max(expression_data) >= 20) {
            seurat <- CreateSeuratObject(counts = as.matrix(expression_data), project = "nn")
            seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
          } else {
            seurat <- CreateSeuratObject(counts = as.matrix(expression_data), project = "nn")
            seurat <- SetAssayData(object = seurat, new.data = as.matrix(expression_data))
          }
          
          seurat <- AddMetaData(object = seurat, metadata = prediction)
          seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
          seurat <- ScaleData(seurat, features = rownames(seurat))
          seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
          seurat <- RunUMAP(seurat, dims = 1:pc_dims, seed.use = seed)
        }
      } else {
        prediction <- copy(cytotrace2_result)
        
        if (max(expression_data) >= 20) {
          seurat <- CreateSeuratObject(counts = as.matrix(expression_data), project = "nn")
          seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
        } else {
          seurat <- CreateSeuratObject(counts = as.matrix(expression_data), project = "nn")
          seurat <- SetAssayData(object = seurat, new.data = as.matrix(expression_data))
        }
        
        seurat <- AddMetaData(object = seurat, metadata = prediction)
        seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
        seurat <- ScaleData(seurat, features = rownames(seurat))
        seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
        seurat <- RunUMAP(seurat, dims = 1:pc_dims, seed.use = seed)
      }
    })
  })
  
  # Merge annotation and predicted labels, if annotation exists
  if (!is.null(annotation)) {
    names(annotation)[1] <- "Phenotype_CT2"
    annotation <- merge(prediction, annotation, by = "row.names", all.x = TRUE)
    rownames(annotation) <- annotation$Row.names
    annotation$Row.names <- NULL
    seurat <- AddMetaData(object = seurat, metadata = annotation["Phenotype_CT2"])
  } else {
    message("Annotation file not provided.")
  }
  
  message("Creating plots.")
  suppressMessages({
    suppressWarnings({
      
      # Plotting
      plot_list <- list()
      labels <- c('Differentiated', 'Unipotent', 'Oligopotent', 'Multipotent', 'Pluripotent', 'Totipotent')
      colors <- c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", "#66C2A5", "#5E4FA2")
      
      x_limits <- range(seurat@reductions$umap@cell.embeddings[,1], na.rm = TRUE)
      y_limits <- range(seurat@reductions$umap@cell.embeddings[,2], na.rm = TRUE)
      
      seurat@meta.data[["CytoTRACE2_Score_clipped"]] <- 5.5 - 6 * seurat@meta.data[["CytoTRACE2_Score"]] 
      seurat@meta.data[["CytoTRACE2_Score_clipped"]] <- -pmax(pmin(seurat@meta.data[["CytoTRACE2_Score_clipped"]], 5), 0)
      
      # UMAP by potency score
      potency_score_umap <- FeaturePlot(seurat, "CytoTRACE2_Score_clipped") +
        scale_colour_gradientn(colours = rev(colors), na.value = "transparent",
                               # breaks=c(0, 0.08333333, 0.25000000, 0.41666667, 0.58333333, 0.75000000, 0.91666667, 1 ),
                               labels = c(labels),
                               limits=c(-5,0), name = "Potency score \n",
                               guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
        xlab("UMAP1") + ylab("UMAP2") + ggtitle("CytoTRACE 2") +
        theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12),
              axis.text = element_text(size = 12), axis.title = element_text(size = 12),   
              plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
              aspect.ratio = 1) +
        coord_cartesian(xlim = x_limits, ylim = y_limits)
      
      plot_list <- c(plot_list, setNames(list(potency_score_umap), "CytoTRACE2_UMAP"))
      
      # UMAP by potency categories
      potency_category_umap <- DimPlot(seurat, reduction = "umap", group.by = "CytoTRACE2_Potency", label = FALSE) +
        scale_color_manual(values = colors, name = "Potency category", breaks = rev(labels)) +
        xlab("UMAP1") + ylab("UMAP2") + ggtitle("CytoTRACE 2") +
        theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12),
              axis.text = element_text(size = 12), axis.title = element_text(size = 12),         
              plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
              aspect.ratio = 1) +
        coord_cartesian(xlim = x_limits, ylim = y_limits)
      
      plot_list <- c(plot_list, setNames(list(potency_category_umap), "CytoTRACE2_Potency_UMAP"))
      
      # UMAP by relative order
      rel_order_umap <- FeaturePlot(seurat, "CytoTRACE2_Relative") +
        scale_colour_gradientn(colours = c("#000004FF", "#3B0F70FF", "#8C2981FF", "#DE4968FF", "#FE9F6DFF", "#FCFDBFFF"),
                               na.value = "transparent", limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
                               labels = c("0.0 (More diff.)", "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"),
                               name = "Relative\norder \n", guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
        xlab("UMAP1") + ylab("UMAP2") + ggtitle("CytoTRACE 2") +
        theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12),
              axis.text = element_text(size = 12), axis.title = element_text(size = 12),     
              plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
              aspect.ratio = 1) +
        coord_cartesian(xlim = x_limits, ylim = y_limits)
      
      plot_list <- c(plot_list, setNames(list(rel_order_umap), "CytoTRACE2_Relative_UMAP"))
      
      # If annotation file is provided
      if (!is.null(annotation)) {
        # UMAP colored by phenotype labels
        phenotype_umap <- DimPlot(seurat, reduction = "umap", group.by = "Phenotype_CT2") +
          xlab("UMAP1") + ylab("UMAP2") + ggtitle("Phenotypes") +
          theme(legend.text = element_text(size = 8), legend.title = element_text(size = 12),
                axis.text = element_text(size = 10), axis.title = element_text(size = 10),            
                plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
                aspect.ratio = 1) +
          coord_cartesian(xlim = x_limits, ylim = y_limits)
        
        plot_list <- c(plot_list, setNames(list(phenotype_umap), "Phenotype_UMAP"))
        
        # Boxplot by phenotype and potency score
        mtd <- seurat@meta.data[c("Phenotype_CT2", "CytoTRACE2_Score")]
        medians <- mtd %>%
          group_by(Phenotype_CT2) %>%
          dplyr::summarise(CytoTRACE2_median_per_pheno = median(CytoTRACE2_Score, na.rm = TRUE)) %>%
          dplyr::arrange(dplyr::desc(CytoTRACE2_median_per_pheno))
        
        mtd <- mtd %>%
          inner_join(medians, by = "Phenotype_CT2")
        
        mtd$Phenotype_CT2 <- factor(mtd$Phenotype_CT2, levels = medians$Phenotype_CT2)
        
        potencyBoxplot_byPheno <- ggplot(mtd[!is.na(mtd$Phenotype_CT2), ], aes(x = Phenotype_CT2, y = CytoTRACE2_Score)) +
          geom_boxplot(aes(fill = CytoTRACE2_median_per_pheno), width = 0.8, alpha = 0.5, outlier.shape = NA) +
          geom_jitter(aes(fill = CytoTRACE2_median_per_pheno), width = 0.05, height = 0, alpha = 0.5, shape = 21, stroke = 0.1, size = 1) +
          theme_classic() +
          scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1),
                             sec.axis = sec_axis(trans = ~., breaks = seq(0, 1, by = 1/12),
                                                 labels = c("", 'Differentiated', "", 'Unipotent', "", 'Oligopotent', "", 'Multipotent', "", 'Pluripotent', "", 'Totipotent', ""))) +
          scale_fill_gradientn(colors = rev(colors), breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1), labels = labels) +
          scale_color_gradientn(colors = rev(colors), breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), limits = c(0, 1), labels = labels) +
          scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
          labs(x = "Phenotype", y = 'Potency score') +
          ggtitle("Developmental potential by phenotype") +
          theme(legend.position = "None", axis.text = element_text(size = 8), axis.text.x = element_text(angle = 90, vjust = 0.5),
                axis.title = element_text(size = 12), legend.text = element_text(size = 12),
                plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
                axis.ticks.y.right = element_line(color = c("black", NA, "black", NA, "black", NA, "black", NA, "black", NA, "black", NA, "black")),
                aspect.ratio = 0.8, axis.ticks.length.y.right = unit(0.3, "cm"))
        
        plot_list <- c(plot_list, setNames(list(potencyBoxplot_byPheno), "CytoTRACE2_Boxplot_byPheno"))
      }
    })
  })
  
  message("Done. You can access any plot directly from the returned list by '$' operator (i.e. plots$CytoTRACE2_Potency_UMAP).")
  return(plot_list)
}