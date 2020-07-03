#' Fetch specified data from TCGA database
#'
#' Fetch specified data from TCGA database by the function 'getFirehoseData' of the package 'RTCGAToolbox'.
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


