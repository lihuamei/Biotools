#' Convert probe IDs to gene IDs
#'
#' Convert a list of probe IDs to specified gene ID type.
#'
#' @param probeIDs A vector of probes.
#' @param targetType Converted targte gene IDs, default: SYMBOL.
#' @return Convert results data.frame.
#' @export probe2GeneIDs
#' @examples
#'
#' convertDf <- probe2GeneIDs(c('1367452_at', '1367453_at', '1367454_at'), targetType = 'SYMBOL')

probe2GeneIDs <- function(probeIDs, targetType = c('SYMBOL')) {
	targetType <- match.arg(targetType)
	convertRes <- switch(
			EXPR    = targetType,
			'SYMBOL' = getSYMBOL(probeIDs, 'hgu133plus2.db'),
		)
	return(convertRes)
}
