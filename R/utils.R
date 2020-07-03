#' Extract subset data frames with the same row or column names by the same order.
#'
#' @param exprSet a expression matrix, which columns are sample,rows are HUGO gene symbols, or probeset ID .
#' @param groupList a vector,as long as the col number for the expression matrix,which describe the group for the samples in exprSet
#' @param prefix The prefix for gct file and cls files.
#' @param destdir where to store the files just download.
#' @return A list of writed 2 files which are the input for GSEA (gct and cls format)
#' * `gctFile` - The saved path of GCT format file
#' * `clsFile` - The saved path of CLS format file
#' @export createGSEAinput
#' @examples
#' 
#' createGSEAinput(exprSet = matrix(rnorm(1000), ncol = 10, nrow = 100), groupList = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2), 'GSE1009')

intersectDfsNames <- function(..., MARGIN = c(1, 2)) {
    dataList <- list(...)
    transBool <- ifelse(match.arg(MARGIN) == 1, FALSE, TRUE)
    
    Reduce(intersect, lapply(dataList, rownames))
}

#' Check numeric data is log2-transformed or not based on the quantile method.
#'
#' @param X A data.frame or a vector of numeric data.
#' @return Bool value. TRUE indicates the checked data is log-transformed, FALSE is not.
#' @export logTransform
#' @examples
#' 
#' logTransform(runif(n = 100, min = 1, max = 1000))

logTransform <- function(X) {
    qx <- as.numeric(quantile(X, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    logBV <- (qx[5] > 100) ||
        (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    return(!logBV)
}

#' Show running information.
#' 
#' @param Message needs to be shown.
#' @param level Specify the output level of the message, 'INFO', 'WARN' or 'ERROR', default: 'INFO'.
#' @param verbose Show running infos or not, default: TRUE.
#' @examples
#'
#' logTransform(runif(n = 100, min = 1, max = 1000))

showRunInfos <- function(msg, level = c('INFO', 'WARN', 'ERROR'), verbose = TRUE) {
	level <- match.arg(level)
	if (verbose) message(sprintf('>> [%s] %s', level, msg))
	if (level == 'ERROR') stop()
}
