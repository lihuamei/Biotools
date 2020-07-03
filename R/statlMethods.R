#' Check numeric data needs log2-transformed or not based on the quantile method.
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
    return(logBV)
}