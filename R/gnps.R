#' Access GNPS information 
#'
#' This function uses GNPS task id and optionally ms2lda task id to retrieve information from ProteoSAFe api.
#' @param gnpsID task id for MZmine workflow.
#' @param ms2ldaID task id for ms2lda workflow.
#' 
#' @return list of dataframes with the results from workflows.
#' 
#' @keywords gnps 
#' @export
#' @examples

access_gnps <- function(gnpsID, ms2ldaID=NULL){
    #url_to_param = paste0('http://proteomics2.ucsd.edu/ProteoSAFe/DownloadResultFile?task=', args[2], '&block=main&file=params/')  
    #nappar <- xmlToList(url_to_param)
    #gtaskid <- unlist(lapply(nappar, function(x) if (x$.attrs=='JOBID') x$text ))    
    base_url <- "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task="
    url_to_metadata = paste0(base_url,  gnpsID, "&block=main&file=metadata_table/")
    url_to_features = paste0(base_url,  gnpsID, "&block=main&file=quantification_table/")
    
    cat('Downloading MZmine feature table...\n')
    feat <- read.csv(url_to_features, check.names=FALSE, stringsAsFactors=FALSE)
    cat('Downloading metadata table...\n')
    meta <- read.delim(url_to_metadata, check.names=FALSE, stringsAsFactors=FALSE)

    if (!is.null(ms2ldaID)) {
        url_to_ms2lda_nodes = paste0(base_url,  ms2ldaID, "&block=main&file=output_results/output_ms2lda_nodes.tsv")
        url_to_motifs_in_scans = paste0(base_url,  ms2ldaID, "&block=main&file=output_results/output_motifs_in_scans.tsv")
        url_to_ms2lda_edges = paste0(base_url,  ms2ldaID, "&block=main&file=output_results/output_ms2lda_edges.tsv")
    
        cat('Downloading ms2lda nodes table...\n')
        ms2lda_nodes <- read.delim(url_to_ms2lda_nodes, check.names=FALSE, stringsAsFactors=FALSE)
        cat('Downloading ms2lda motifs table...\n')
        motifs_in_scans <- read.delim(url_to_motifs_in_scans, check.names=FALSE, stringsAsFactors=FALSE)
        cat('Downloading ms2lda edges table...\n')
        ms2lda_edges <- read.delim(url_to_ms2lda_edges, check.names=FALSE, stringsAsFactors=FALSE)
    } else {
        ms2lda_nodes <- NULL 
        motifs_in_scans <- NULL 
        ms2lda_edges <- NULL 
    }
    return(list(features=feat, metadata=meta, ms2lda_nodes=ms2lda_nodes, motifs_in_scans=motifs_in_scans, ms2lda_edges=ms2lda_edges))
}
