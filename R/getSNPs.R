#' Get list of SNPs from EBI's GWAS server
#'
#' @param queryTerm Query term 
#'
#' @return Returns a data fram of SNPs from EBI's server
#' @export
#' @import httr
#' @importFrom utils URLdecode URLencode
#'
#' @examples
#' getSNPs('breast cancer')


getSNPs <- function(queryTerm) {
  library(httr)
  url <- parse_url('https://www.ebi.ac.uk/gwas/api/search/downloads?q=&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&traitfilter%5B%5D=&dateaddedfilter=&facet=association')
  url$query$q <- URLdecode(paste0('text:%22',URLencode(queryTerm),'%22'))
  url <- build_url(url)
  dat <- content(GET(url), as = 'text', encoding = 'UTF-8')
  tmp <- unlist( strsplit( dat, "\r\n", fixed = TRUE ) )
  ## This following split_to_df function is originally from rsnps package but modified to my needs.
  split_to_df <- function(x, sep, fixed=FALSE, perl=TRUE, useBytes=FALSE, names=NULL) {
    
    x <- as.character(x)
    
    if( fixed ) {
      perl <- FALSE
    }
    
    tmp <- strsplit( x, sep, fixed=fixed, perl=perl, useBytes=useBytes )
    if( length( unique( unlist( lapply( tmp, length ) ) ) ) > 1 ) {
      stop("non-equal lengths for each entry of x post-splitting")
    }
    tmp <- unlist( tmp )
    tmp <- as.data.frame( 
      matrix( tmp, ncol = (length(tmp) / length(x)), byrow=TRUE ),
      stringsAsFactors=FALSE, optional=TRUE 
    )
    
    if( !is.null(names) ) {
      names(tmp) <- names
    } else {
      names(tmp) <- tmp[1,]
      tmp <- tmp[-1,]
    }
    
    return(tmp)
  }
  out <- split_to_df( tmp, sep = "\t", fixed = TRUE )
}
