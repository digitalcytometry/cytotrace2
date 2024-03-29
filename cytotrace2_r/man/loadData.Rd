% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{loadData}
\alias{loadData}
\title{Load Data from File}
\usage{
loadData(filepath)
}
\arguments{
\item{filepath}{The path to the file containing the data
(count matrix, samples as columns, genes as rows)}
}
\value{
The loaded data as a data frame.
}
\description{
This function loads data from a specified file path.
}
\examples{
data <- loadData("/path/to/file.tsv")

}
