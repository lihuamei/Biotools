#' Fetch genomic datasets based on the given cancer types.
#'
#' Fetch the genomic datasets for the given cancer type from TCGA database using the package 'TCGAbiolinks'. All cancer type names can be obtained from https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html.
#'
#' @param cancerTypes Specify the cancer types for downloading the genimic datasets from the TCGA database.
#' @param keep Keep downloaded files and directories or not, default: FALSE.
#' @param ... Other parameters can be obtained through GDCquery functions.
#' @return A list of download datasets for the given cancer types.
#' * `data` - Fetched data matrix.
#' * `clinical` - Clinical information corressponding to the cancer type.
#' @export downloadTcgaData
#' @examples
#' 
#' results <- downloadTcgaData(cancerTypes = c('TCGA-GBM'), data.category = 'Gene expression', data.type = 'Gene expression quantification', platform = 'Illumina HiSeq', file.type  = 'normalized_results', experimental.strategy = 'RNA-Seq', legacy = TRUE)

downloadTcgaData <- function(cancerTypes = getGDCprojects()$project_id[grepl('TCGA', getGDCprojects()$project_id)], keep = FALSE, ...) {
	lapply(cancerTypes, function(cancer) {
		query <- GDCquery(project = cancer, ...)
		GDCdownload(query)
		dataTmp <- GDCprepare(query)
		if (!keep) unlink(cancerTypes, recursive = TRUE)
		return(list(data = assay(dataTmp), clinical = colData(dataTmp)))
	}) -> downloadedResults
	names(downloadedResults) <- cancerTypes
	if (!keep) {
		file.remove('MANIFEST.txt')
		unlink('GDCdata', recursive = TRUE)
	}
	return(downloadedResults)
}
