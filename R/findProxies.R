#' Submit list of SNPs to SNPsnap server
#'
#' @param path2RefDir Path to reference directory, must contain the prefix of the plink reference files bed, bim and fam
#' @param path2leadSNPList Full path to the list of lead snps
#' @param proxyWindow --ld-window-kb;  default is 1000
#' @param proxyCorr plink parameter --ld-window-r2 ; default is 0.8
#' @param path2Proxies patht to directory where the SNPs are
#' @return
#' @export
#' @examples
#' findProxies()

findProxies <- function(path2RefDir,path2leadSNPList,proxyWindow=1000,proxyCorr=0.8,path2Proxies)
  {
  if (! is.character(path2RefDir))  {stop("Parameter path2RefDir should be a string.", call. = FALSE)}
  if ( !all( sapply(paste(path2RefDir,c('bed','bim','fam'),sep  = '.'),function(x) (file.exists(x))))) {
    stop("contain the prefix of the plink reference files bed, bim and fam", call. = FALSE) 
  }
  if (! is.character(path2leadSNPList))  {stop("Parameter path2RefDir should be a string.", call. = FALSE)}
  if (! file.exists(path2leadSNPList))  {stop("Lead snp file does not exist.", call. = FALSE)}
  if (! is.character(path2Proxies))  {stop("Parameter path2Proxies should be a string.", call. = FALSE)}
  if (! dir.exists(path2Proxies)) { dir.create(path2Proxies, recursive = T)}
  if (proxyWindow%%1!=0 || proxyWindow < 1) {stop("Parameter proxyWindow should be positive integer,", call. = FALSE)}
  if (proxyCorr < -1.0 || proxyCorr > 1) {stop("Invalid correlation.", call. = FALSE)}
  commandstr <- paste('plink','--bfile',  path2RefDir,
                              '--ld-snp-list', path2leadSNPList,
                              '--ld-window-kb', as.character(proxyWindow),
                              '--ld-window-r2',as.character(proxyCorr),
                              '--r2',
                              '--out',as.character(path2Proxies))
  
  system(commandstr)
  
}