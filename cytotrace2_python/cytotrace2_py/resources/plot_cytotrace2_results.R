
suppressWarnings({suppressPackageStartupMessages(library(argparse))})
suppressWarnings({suppressPackageStartupMessages(library(Seurat))})
suppressWarnings({suppressPackageStartupMessages(library(ggplot2))})
suppressWarnings({suppressPackageStartupMessages(library(stringr))})
suppressWarnings({suppressPackageStartupMessages(library(dplyr))})


loadData <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("File not found.")
  }

  # load data
  data <- data.table::fread(filepath, header = TRUE, sep = "\t", check.names = FALSE, data.table = FALSE)

  # check for uniqueness of gene names
  if (any(duplicated(data[[1]]))) {
    stop("Please make sure the gene names are unique")
  }
  
  row.names(data) <- data[[1]]
  row.names(data) <- gsub(".", "-", row.names(data), fixed = T)
  data[[1]] <- NULL
  
  # check for uniqueness of gene names
  if (any(duplicated(colnames(data)))) {
    stop("Please make sure the cell/sample names are unique in the expression data.")
  }

  return(data)
}


plotData <- function(expression_path,
                     result_path,
                     annotation_path,
                     plot_dir,
                     num_pcs = 30,
                     seed = 14) {

  set.seed(seed)
  
  expression_data <- loadData(expression_path)
  prediction <- loadData(result_path)
  expression_data <- expression_data[,rownames(prediction)]

  dir.create(plot_dir,showWarnings=F,recursive=T)

  if(annotation_path != ""){
    annotation <-  data.table::fread(annotation_path, header = TRUE, sep = "\t", check.names = FALSE, data.table = FALSE)
    if (any(duplicated(annotation[[1]]))) {
      stop("Please make sure the cell/sample names are unique in annotation data.")
    }

    row.names(annotation) <- annotation[[1]]
    annotation[[1]] <- NULL

  } else{
    annotation <- NULL
  }


  message("Preparing input for visualization.")
  suppressMessages({
    suppressWarnings({
  
    # Normalizing data - LogNormalize
    # check if not already log-normalized
    # rule-of-thumb max value of the matrix should not be less than 20
    if (!(max(expression_data) < 20)) {
      seurat <- CreateSeuratObject(counts = as.matrix(expression_data), project = "nn", min.cells = 3, min.features = 200)
      seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000) ##defaults
    } else {
      seurat <- CreateSeuratObject(counts = as.matrix(expression_data),  project = "nn", min.cells = 3, min.features = 200)
      seurat <- SetAssayData(object = seurat, new.data = as.matrix(expression_data))
    }
 

    invisible(capture.output(
    # Calculate variable features
    seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
    )) 
      
    # Scale data
    all.genes <- rownames(seurat)
    seurat <- ScaleData(seurat, features = all.genes)

    # PCA
    seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))

    invisible(
      capture.output(
        # visualize clustering
        seurat <- RunUMAP(seurat, dims = 1:num_pcs) # default to 30
      )
    )
    }) })

  # Add needed values to metadata
  # merge annotation and predicted labels, if annotation exists
  # if annotation file is provided
  if (annotation_path != "")  {

    names(annotation)[1] <- "Phenotype"
    annotation <- merge(prediction, annotation, by = "row.names", all.x = TRUE)
    rownames(annotation) <- annotation$Row.names
    annotation$Row.names <- NULL
    
    seurat <- AddMetaData(
      object = seurat,
      metadata = prediction
    )
    
    seurat <- AddMetaData(
      object = seurat,
      metadata = annotation["Phenotype"]
    )
    
  } else {
    
    seurat <- AddMetaData(
      object = seurat,
      metadata = prediction
    )
  } 


  suppressMessages({
    suppressWarnings({
  # plotting
  labels <- c('Differentiated', 'Unipotent', 'Oligopotent', 'Multipotent', 'Pluripotent', 'Totipotent')
  colors <- c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", "#66C2A5", "#5E4FA2")

  # axis limits
  x_limits <- range(seurat@reductions$umap@cell.embeddings[,1], na.rm = TRUE)
  y_limits <- range(seurat@reductions$umap@cell.embeddings[,2], na.rm = TRUE)

  # UMAP by potency score
  seurat@meta.data[["CytoTRACE2_Score_clipped"]] <- 5.5 - 6 * seurat@meta.data[["CytoTRACE2_Score"]] 
  
   # Then, we apply the minimum and maximum cutoffs
  seurat@meta.data[["CytoTRACE2_Score_clipped"]]  <- -pmax(pmin(seurat@meta.data[["CytoTRACE2_Score_clipped"]], 5), 0)
  
  # UMAP by potency score
  potency_score_umap <- FeaturePlot(seurat, "CytoTRACE2_Score_clipped") +
          scale_colour_gradientn(colours = rev(colors), na.value = "transparent",
                                 # breaks=c(0, 0.08333333, 0.25000000, 0.41666667, 0.58333333, 0.75000000, 0.91666667, 1 ),
                                 labels = c(labels),
                                 limits=c(-5,0), name = "Potency score \n",
                                 guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    ggtitle("CytoTRACE 2") +
    theme(legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),   
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20))) +
    theme(aspect.ratio = 1)  +
    coord_cartesian(xlim = x_limits, ylim = y_limits)


  pdf(paste0(plot_dir,'/CytoTRACE2_Score_UMAP.pdf'))
  print(potency_score_umap)
  dev.off()

  # UMAP by potency categories
  potency_category_umap <- DimPlot(seurat, reduction = "umap", group.by = "CytoTRACE2_Potency", label = FALSE) +
    scale_color_manual(values = colors, name = "Potency category",
                       breaks = rev(c('Differentiated', 'Unipotent', 'Oligopotent', 'Multipotent', 'Pluripotent', 'Totipotent'))) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    ggtitle("CytoTRACE 2") +
    theme(legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          axis.text = element_text(size = 12),
         axis.title = element_text(size = 12),         
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20))) +
    theme(aspect.ratio = 1) +
    coord_cartesian(xlim = x_limits, ylim = y_limits)


  pdf(paste0(plot_dir,'/CytoTRACE2_Potency_UMAP.pdf'))
  print(potency_category_umap)
  dev.off()

  # UMAP by relative order
  rel_order_umap <- FeaturePlot(seurat, "CytoTRACE2_Relative") +
    scale_colour_gradientn(colours = (c( "#000004FF", "#3B0F70FF", "#8C2981FF", "#DE4968FF", "#FE9F6DFF", "#FCFDBFFF")),
                           na.value = "transparent",
                           limits=c(0,1),
                           breaks = seq(0,1, by = 0.2),
                           labels=c("0.0 (More diff.)", "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"),
                           name = "Relative\norder \n" ,
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    ggtitle("CytoTRACE 2") +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme(legend.text = element_text(size = 10),
          legend.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),     
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20))) +
    theme(aspect.ratio = 1) +
    coord_cartesian(xlim = x_limits, ylim = y_limits)

  pdf(paste0(plot_dir,'/CytoTRACE2_Relative_UMAP.pdf'))
  print(rel_order_umap)
  dev.off()

  # if annotation file is provided
  if (annotation_path != "")  {
    # UMAP colored by phenotype labels
    phenotype_umap <- DimPlot(seurat, reduction = "umap", group.by = "Phenotype", label = FALSE) +
      xlab("UMAP1") +
      ylab("UMAP2") +
      ggtitle("Phenotypes") +
      theme(legend.text = element_text(size = 8),
            legend.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 10),            
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20))) +
      theme(aspect.ratio = 1) +
      coord_cartesian(xlim = x_limits, ylim = y_limits)

    pdf(paste0(plot_dir,'/CytoTRACE2_phenotype_UMAP.pdf'))
    print(phenotype_umap)
    dev.off()

    # Boxplot by phenotype and potency_score
    # adjust ordering of boxes according to the potency score range
    mtd <- seurat@meta.data[c("Phenotype", "CytoTRACE2_Score")]
    medians <- mtd %>%
      group_by(Phenotype) %>%
      summarise(CytoTRACE2_median_per_pheno = median(CytoTRACE2_Score, na.rm = TRUE)) %>%
      arrange(desc(CytoTRACE2_median_per_pheno))
    
    # Join the median scores back to the original dataframe for coloring
    mtd <- mtd %>%
      inner_join(medians, by = "Phenotype")
    
    # Order Phenotype by median CytoTRACE2_Score for plotting
    mtd$Phenotype <- factor(mtd$Phenotype, levels = medians$Phenotype)
    
    
    # Create a boxplot for CytoTRACE2_Score by Phenotype

    potencyBoxplot_byPheno <- ggplot(mtd[!is.na(mtd$Phenotype), ], aes(x = Phenotype, y = CytoTRACE2_Score)) +
      geom_boxplot(aes(fill = CytoTRACE2_median_per_pheno), width = 0.8, alpha = 0.5, outlier.shape = NA) +
      geom_jitter(aes(fill = CytoTRACE2_median_per_pheno), width = 0.05, height = 0, alpha = 0.5, shape = 21, stroke = 0.1, size = 1) +
      theme_classic() +
      scale_y_continuous(breaks=seq(0, 1, by = 0.2),  limits = c(0,1),
                         sec.axis = sec_axis(trans = ~., breaks = seq(0, 1, by = 1/12),
                                             labels = c("", 'Differentiated', "",'Unipotent', "",'Oligopotent', "",'Multipotent',"", 'Pluripotent', "",'Totipotent', "") )) +
      scale_fill_gradientn(colors = rev(colors),
                           breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                           limits = c(0,1), 
                           labels = c(labels))+
      scale_color_gradientn(colors = rev(colors),
                            breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                            limits = c(0,1), 
                            labels = c(labels))+
      scale_x_discrete(labels = function(x) str_wrap(x, width = 10))  +
      labs(x = "Phenotype", y = 'Potency score') +
      ggtitle("Developmental potential by phenotype") +
      theme(legend.position = "None",
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 12),         
            legend.text = element_text(size = 12),
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
            axis.ticks.y.right  = element_line(color = c("black", NA, "black",NA, "black",NA, "black",NA,"black", NA, "black",NA, "black")),
            aspect.ratio = 0.8,
            axis.ticks.length.y.right  = unit(0.3, "cm"))

    pdf(paste0(plot_dir,'/CytoTRACE2_box_plot_by_phenotype.pdf'))
    print(potencyBoxplot_byPheno)
    dev.off()
  } }) })
  
}


parser <- ArgumentParser()
parser$add_argument("--expression-path", type="character", help="Path to file containing input expression data")
parser$add_argument("--result-path", type="character", help="Path to file containing CytoTRACE2 results")
parser$add_argument("--annotation-path", type="character", default="", help="Path to file containing input data annotation")
parser$add_argument("--plot-dir", type="character", default="", help="Directory to which to save plots")
parser$add_argument("--num-pcs", type="integer", default=30, help="Integer, indicating the number of principal components to use for Seurat PCA")
parser$add_argument("--seed", type="integer", default=14, help="Integer, specifying the seed for reproducibility in random processes")

args <- parser$parse_args()
plotData(expression_path = args$expression_path, result_path = args$result_path, annotation_path = args$annotation_path, plot_dir = args$plot_dir, num_pcs = args$num_pcs, seed = args$seed)
