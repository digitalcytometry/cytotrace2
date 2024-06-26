% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postprocessing.R
\name{smoothDatakNN}
\alias{smoothDatakNN}
\title{Smooth Data with kNN}
\usage{
smoothDatakNN(ranked_df, predicted_df, top_genes, max_pcs, seed)
}
\arguments{
\item{ranked_df}{dataframe of ranked gene expression (ranked data, transposed) where rows are genes and columns are cells}

\item{predicted_df}{dataframe of 4 columns containing CytoTRACE 2 model predicted and postprocessed potency scores and categories}

\item{top_genes}{list of the top 1000 most variable genes (output of dispersion function)}

\item{max_pcs}{integer indicating the maximum number of principal components to use in the smoothing by kNN step (default is 200).}

\item{seed}{integer specifying the seed for reproducibility in random processes (default is 14).}
}
\value{
A dataframe with final updated predictions
}
\description{
This applies k-nearest neighbor smoothing to datasets with >100 cells.
It identifies the nearest neighbors of each cell based on top principal components,
and assigns the median score of the cell and its  nearest neighbors, as the final
potency score for the cell, and updates potency category accordingly.
}
\examples{
 ranked_data --> first element in the output of preprocessData()
 predicted_df --> output of binData()
 ranked_df <- data.frame(t(ranked_data))
 rownames(ranked_df) <- row.names(expression_data)
 colnames(ranked_df) <- colnames(expression_data)
 predicted_df <- smoothDatakNN(ranked_df, predicted_df, top_genes)
}
