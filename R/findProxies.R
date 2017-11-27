#' Submit list of SNPs to SNPsnap server
#'
#' @param path2PlinkPrefix Path to reference directory, must contain the prefix of the plink reference files bed, bim and fam
#' @param path2leadSNPList Full path to the list of lead snps
#' @param proxyWindow --ld-window-kb;  default is 1000
#' @param proxyCorr plink parameter --ld-window-r2 ; default is 0.8
#' @param path2Proxies patht to directory where the SNPs are
#' @return
#' @export
#' @examples
#' findProxies()

findProxies <- function(path2PlinkPrefix,snplist,proxyWindow=1000,proxyCorr=0.8,path2Proxies)
  {
  if (! is.character(path2PlinkPrefix))  {stop("Parameter path2PlinkPrefix should be a string.", call. = FALSE)}
  if ( !all( sapply(paste(path2PlinkPrefix,c('bed','bim','fam'),sep  = '.'),function(x) (file.exists(x))))) {
    stop("Should contain the prefix of the plink reference files bed, bim and fam", call. = FALSE) 
  }
  if (! is.vector(snplist))  {stop("Parameter path2RefDir should be a vector.", call. = FALSE)}
  if (! is.character(path2Proxies))  {stop("Parameter path2Proxies should be a string.", call. = FALSE)}
  if (! dir.exists(dirname(path2Proxies))) { dir.create(dirname(path2Proxies), recursive = T)}
  if (proxyWindow%%1!=0 || proxyWindow < 1) {stop("Parameter proxyWindow should be positive integer,", call. = FALSE)}
  if (proxyCorr < -1.0 || proxyCorr > 1) {stop("Invalid correlation.", call. = FALSE)}
  
  tempfilename=tempfile()
  write.table(as.data.frame(unique(snplist)),file=tempfilename,quote=F,sep="\t",row.names=F,col.names=F)
  
  commandstr <- paste('plink','--bfile',  path2PlinkPrefix,
                              '--ld-snp-list', tempfilename,
                              '--ld-window-kb', as.character(proxyWindow),
                              '--ld-window-r2',as.character(proxyCorr),
                              '--r2',
                              '--out',as.character(path2Proxies))
  
  system(commandstr)
  
}