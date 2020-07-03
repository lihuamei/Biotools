#' Download GEO datasets from Gene Expression Omnibus databse
#'
#' Download a list of GEO datasets from Gene Expression Omnibus databse.
#'
#' @param geoIDs A vector of GSE accession number.
#' @param verbose Show running infos or not, default: TRUE.
#' @param log2Trans log2(X + 1) transfer for the expression data, default: FALSE.
#' @return A list of downloaded gene expression profiles of each GSE accession number.
#' @export downGEODatasets
#' @examples
#' 
#' exprs.lst <- downGEODatasets(geoIDs = c('GSE19830'), log2Trans = FALSE)

downGEODatasets <- function(geoIDs, log2Trans = FALSE, verbose = TRUE) {
	sapply(geoIDs, function(gse) {
		showRunInfos(msg = sprintf('Fetching %s', gse), level = 'INFO', verbose = verbose)
		status <- try(gsets <- getGEO(gse, GSEMatrix = TRUE, AnnotGPL = TRUE))
		
		if (class(status) == 'try-error') {
			showRunInfos(msg = sprintf('Fetching %s failed', gse), level = 'WARN', ...)
			next
		}
		lapply(1 : length(gsets), function(idx){
			gset  <- gsets[[idx]]
			fvarLabels(gset) <- make.names(fvarLabels(gset))
			ex <- exprs(gset)
			logBool <- logTransform(ex)
			if (log2Trans && !logBool) { ex[which(ex <= 0)] <- NaN; ex <- log2(1 + ex) }
			return(ex)
		}) -> exprs
	}) -> fetched.exprs
	return(fetched.exprs)
}
