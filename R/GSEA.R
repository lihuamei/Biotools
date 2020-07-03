#' create gct file and cls files
#'
#' create gct file and cls file for GSEA according the expression matrix and group information
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

createGSEAinput <- function(exprSet, groupList, prefix = 'GSEA', destDir = '.') {
    gctFile <- file.path(destDir, paste0(prefix, '.gct'))
    sink(gctFile)
    cat('#1.2\n')
    cat(paste0(nrow(exprSet), "\t", length(groupList), "\n"))
    sink()
    gctOut <- cbind(symbol = rownames(exprSet), description = "na", exprSet)
    write.table(gctOut, gctFile, append = TRUE, quote = FALSE, row.names = FALSE, sep = '\t')

    clsFile = file.path(destDir, paste0(prefix, '.cls'))
    sink(clsFile)
    cat(paste0(length(groupList), " ", length(unique(groupList)), " 1\n"))
    cat(paste0("# ", paste(unique(groupList), collapse = " "), "\n"))
    cat(paste(groupList, collapse = " "))
    sink()
    return(list(gctFile = gctFile, clsFile = clsFile))
}

#' Hypergeometric Tests for GO/KEGG test
#'
#' Given GeneID2Path,Path2GeneID, diff_gene,universeGeneIds, this function will compute Hypergeometric Tests for each path (GO/KEGG),
#' It can't hold on the structure of the GO graph, just a simple path.
#'
#' @param diffExprGenes   a vector which contain the significantly DEG list.
#' @param geneID2Path a list which one entrez gene id to multiple pathway id.
#' @param path2GeneID a list which one pathway id to multiple entrezgene id.
#' @param universeGeneIds a vector which contain the backgroud gene list, probably 20,000 genes.
#' @return a data.frame, each row is a Hypergeometric Tests result for each pathway.
#' @export hyperGtestTerms
#' @examples
#' 
#' 

hyperGtestTerms <- function(diffExprGenes, geneID2Path, path2GeneID, universeGeneIds) {
    diffExprGenes <- unique(diffExprGenes); universeGeneIds <- unique(universeGeneIds)
    diffGeneHasPathes <- intersect(diffExprGenes, names(geneID2Path))
    n <- length(diff_gene)
    N <- length(universeGeneIds)
    
    results <- sapply(X = names(Path2GeneID), FUN = function(genes) {
        M <- length(intersect(genes, universeGeneIds))
        if (M < 5) next      
        k <- intersect(diffGeneHasPathes, genes) %>% length
        if (k == 0) next
        expCount <- n * M / N
        OddsRatio <- k / expCount
        pVal <- phyper(k - 1, M, N - M, n, lower.tail = F)
        return(c(i, p, OddsRatio, expCount, k, M))
    }) %>% data.frame(., stringsAsFactors = FALSE)

    colnames(results) <- c('PathwayID', 'Pvalue', 'OddsRatio', 'ExpCount', 'Count', 'Size')
    results$p.adjust <- p.adjust(p = results$Pvalue, method = 'BH')
    results <- results[order(results$Pvalue), ]
    rownames(results) <- 1 : nrow(results)
    return(results)
}