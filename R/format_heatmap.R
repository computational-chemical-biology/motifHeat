#' Wrapper to format heamap from defined GNPS-MZmine inputs 
#'
#' This function allows you to subset tables, create multiple colored bars for samples, change color scales, etc.
#' @param featureTable MZmine format feature table.
#' @param metadataTable qiime2 format metadata table.
#' @param selectField metadata field to subselect. Default is NULL.
#' @param selectValue metadata field value to subselect. Default is NULL.
#' @param factorColList nested list containing lists of: a vector of colors, one for each factor level, and a string with a metadata field name. Default is list().
#' @param  colorScale hexadecimal vector of colors. Default is heat.colors(75).
#' @param  main title for heatmap plot. Default is ''.
#' 
#' @return  heatmap list.
#' 
#' @keywords heatmap 
#' @export
#' @examples
#' @import gplots

format_heatmap <- function(featureTable, metadataTable, selectField=NULL, selectValue=NULL,
			   factorColList=list(), colorScale=rev(heat.colors(75)), main='', ...){ 
    if (!is.null(selectField)) {
        meta2 <- meta[meta[,selectField]==selectValue,] 
    } else {
        meta2 <- meta
    }
    
    # Format MZmine feature table
    tab2 <- tab[, grep('Peak area', colnames(tab))] 
    tab2 <- t(tab2) 
    rownames(tab2) <- sub(' Peak area$', '', rownames(tab2)) 
    # order samples according to metadata
    tab2 <- tab2[meta2[,1],] 
    tab2 <- tab2[,apply(tab2, 2, sum)!=0]
   
    colList <- list()
    if (length(factorColList)) { 
        fnames <- c()
	for(i in 1:length(factorColList)) {
            colList[[i]] <- factorColList[[i]]$colors 
            names(colList[[i]]) <- unique(meta2[, factorColList[[i]]$factor])
            colList[[i]] <- colList[[i]][meta2[, factorColList[[i]]$factor]]
            fnames <- append(fnames, factorColList[[i]]$factor)
        }
        rlab <- do.call(rbind, colList)
        rownames(rlab) <- fnames 
    } 
    
    mydist <- function(c) { dist(c, method="canberra") } 
    myclust <- function(c) { hclust(c, method="ward.D") }
    
    h <- heatmap.3(tab2, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="column", dendrogram="both", margins=c(6,12),
     Rowv=TRUE, Colv=TRUE,  RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE,
     density.info="none", trace="none", main=main, labCol=FALSE, labRow=rownames(tab2), cexRow=1, col=colorScale,
     ColSideColorsSize=7, RowSideColorsSize=2, KeyValueName="Ion intensity", ...)
     return(h)
}


