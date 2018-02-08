#' Get list of SNPs from EBI's GWAS server
#'
#' GetSNPs takes a character string and sends it to NHGRI-EBI GWAS Catalog as the query term 
#' and returns a data table containing all data matching to the query term.
#'
#' @param queryTerm The query term for EBI GWAS. See https://www.ebi.ac.uk/gwas/ for details.
#'
#' @return Returns a data frame of SNPs from EBI's server
#' @export
#' @importFrom httr parse_url build_url content GET
#' @importFrom utils URLdecode URLencode
#'
#' @examples
#' getSNPs('breast cancer')


getSNPs <- function(queryTerm) {
  url <- parse_url(gwasURL)
  url$query$q <- URLdecode(paste0('text:%22',URLencode(queryTerm),'%22'))
  url <- build_url(url)
  dat <- content(GET(url), as = 'text', encoding = 'UTF-8')
  tmp <- unlist( strsplit( dat, "\r\n", fixed = TRUE ) )
  out <- split_to_df( tmp, sep = "\t", fixed = TRUE )
  return(out)
}
