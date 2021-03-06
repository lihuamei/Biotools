% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{createPhenotypeClass}
\alias{createPhenotypeClass}
\title{Create phenotype class matrix.}
\usage{
createPhenotypeClass(snList)
}
\arguments{
\item{snList}{A list of sample names.}
}
\value{
Phenotype class matrix, which rows represent the names of snList, columns represent sample names. When sample belongs to class A, the entry = 1, ortherwise, entry = 0.
}
\description{
Create phenotype class matrix.
}
\examples{

createPhenotypeClass(list(A = c('S1', 'S2'), B = c('S3', 'S4'))) 
}
