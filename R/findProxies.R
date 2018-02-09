#' Find proxies to SNPs using plink
#' 
#' ClumpSNPs is a wrapper to plink 1.9. The function computes proxy SNPs for the input snps 
#' and writes them to a file.
#' 
#' @param path2PlinkPrefix Plink bfile parameter. It is the path to reference directory including prefix of the plink reference files bed, bim and fam
#' @param path2leadSNPList Full path to the SNP list 
#' @param ld_window_kb Plink parameter ld_window_kb denotes the maximum distance between LD buddies (default is 1000).
#' @param ld_window_r2 Plink parameter ld_window_r2 denoted the minimum correlation  between LD buddies (default is 0.8).
#' @param path2Proxies Path to directory where the result will be written.
#'  
#' @author Kari Nousiainen 
#' @export
#' 
#' @examples
#' findProxies(path2PlinkPrefix, snplist, path2Proxies)

findProxies <- function(path2PlinkPrefix,snplist,ld_window_kb=1000,ld_window_r2=0.8,path2Proxies)
  {
  if (! is.character(path2PlinkPrefix))  {stop("Parameter path2PlinkPrefix should be a string.", call. = FALSE)}
  if ( !all( sapply(paste(path2PlinkPrefix,c('bed','bim','fam'),sep  = '.'),function(x) (file.exists(x))))) {
    stop("Should contain the prefix of the plink reference files bed, bim and fam", call. = FALSE) 
  }
  if (! is.vector(snplist))  {stop("Parameter path2RefDir should be a vector.", call. = FALSE)}
  if (! is.character(path2Proxies))  {stop("Parameter path2Proxies should be a string.", call. = FALSE)}
  if (! dir.exists(dirname(path2Proxies))) { dir.create(dirname(path2Proxies), recursive = T)}
  if (ld_window_kb%%1!=0 || ld_window_kb < 1) {stop("Parameter ld_window_kb should be positive integer,", call. = FALSE)}
  if (ld_window_r2 < -1.0 || ld_window_r2 > 1) {stop("Invalid correlation.", call. = FALSE)}
  
  tempfilename=tempfile()
  write.table(as.data.frame(unique(snplist)),file=tempfilename,quote=F,sep="\t",row.names=F,col.names=F)
  
  commandstr <- paste('plink','--bfile',  path2PlinkPrefix,
                              '--ld-snp-list', tempfilename,
                              '--ld-window-kb', as.character(ld_window_kb),
                              '--ld-window-r2',as.character(ld_window_r2),
                              '--ld-window 99999',
                              '--r2',
                              '--out',as.character(path2Proxies))
  
  system(commandstr)
  
}