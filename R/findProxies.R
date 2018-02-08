#' Find proxies to SNPs using plink
#' 
#' Computes proxies for given SNPs using plink
#' 
#' @param path2PlinkPrefix Path to reference directory, must contain the prefix of the plink reference files bed, bim and fam
#' @param path2leadSNPList Full path to the list of lead snps
#' @param ld-window-kb --ld-window-kb;  default is 1000
#' @param ld-window-r2 plink parameter --ld-window-r2 ; default is 0.8
#' @param path2Proxies path to directory where the SNPs are
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
  if (ld-window-kb%%1!=0 || ld-window-kb < 1) {stop("Parameter ld-window-kb should be positive integer,", call. = FALSE)}
  if (ld-window-r2 < -1.0 || ld-window-r2 > 1) {stop("Invalid correlation.", call. = FALSE)}
  
  tempfilename=tempfile()
  write.table(as.data.frame(unique(snplist)),file=tempfilename,quote=F,sep="\t",row.names=F,col.names=F)
  
  commandstr <- paste('plink','--bfile',  path2PlinkPrefix,
                              '--ld-snp-list', tempfilename,
                              '--ld-window-kb', as.character(ld-window-kb),
                              '--ld-window-r2',as.character(ld-window-r2),
                              '--ld-window 99999',
                              '--r2',
                              '--out',as.character(path2Proxies))
  
  system(commandstr)
  
}