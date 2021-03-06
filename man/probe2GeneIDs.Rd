% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convertGenes.R
\name{probe2GeneIDs}
\alias{probe2GeneIDs}
\title{Convert probe IDs to gene IDs}
\usage{
probe2GeneIDs(probeIDs, targetType = c("SYMBOL"))
}
\arguments{
\item{probeIDs}{A vector of probes.}

\item{targetType}{Converted targte gene IDs, default: SYMBOL.}
}
\value{
Convert results data.frame.
}
\description{
Convert a list of probe IDs to specified gene ID type.
}
\examples{

convertDf <- probe2GeneIDs(c('1367452_at', '1367453_at', '1367454_at'), targetType = 'SYMBOL')
}
