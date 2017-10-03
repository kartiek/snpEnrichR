#' Submit list of SNPs to SNPsnap server
#'
#' @param path2RefDir Path to reference directory, must contain the prefix of the plink reference files bed, bim and fam
#' @param path2leadSNPList Full path to the list of lead snps
#' @param proxyWindow --ld-window-kb;  default is 1000
#' @param proxyCorr plink parameter --ld-window-r2 ; default is 0.8
#' @param path2Proxies
#' @param useMaf, a logical value to decide whether maf pameter is used;  default is TRUE
#' @param maf plink parameter maf;  default is 0.01
#' @param useHWE, a logical value to decide whether HWE pameter is used;  default is TRUE
#' @param hwe plink parameter hwe;  default is 10e-6
#' @param useGENO, a logical value to decide whether geno pameter is used;  default is TRUE
#' @param geno plink parameter 0.1; default is  0.1 
#' @return Returns a URL from which results can be downloaded
#' @export
#' @examples
#' findProxies()

findProxies <- function(path2RefDir,ldSNPList,proxyWindow=1000,proxyCorr=0.8,path2Proxies,useMaf=TRUE,maf=0.01,useHWE=TRUE,hwe=10e-6,useGENO=TRUE,geno=0.1)
  {
  if (! is.character(path2RefDir))  {stop("Parameter path2RefDir should be a string.", call. = FALSE)}
  if ( !all( sapply(paste(path2RefDir,c('bed','bim','fam'),sep  = '.'),function(x) (file.exists(x))))) {
    stop("contain the prefix of the plink reference files bed, bim and fam", call. = FALSE) 
  }
  if (! is.character(path2leadSNPList))  {stop("Parameter path2RefDir should be a string.", call. = FALSE)}
  if (! file.exists(path2leadSNPList))  {stop("Lead snp file does not exist.", call. = FALSE)}
  if (! is.character(path2Proxies))  {stop("Parameter path2Proxies should be a string.", call. = FALSE)}
  if ( dir.exists(path2Proxies)) { dir.create(path2Proxies)}
  if (proxyWindow%%1!=0 || proxyWindow < 1) {stop("Parameter proxyWindow should be positive integer,", call. = FALSE)}
  if (proxyCorr < -1.0 || proxyCorr > 1) {stop("Invalid correlation.", call. = FALSE)}
  if (! is.logical(useMaf) ) {stop("Parameter useMaf must be logical.", call. = FALSE)}
  if (! is.numeric(maf)) {stop("Parameter maf must be numeric.",call=FALSE)} 
  if (! is.logical(useHWE) ) {stop("Parameter useHWE must be logical.", call. = FALSE)}
  if (! is.numeric(hwe) ) {stop("hwe must be numeric.",call=FALSE)}
  if (! is.logical(useGENO) ) {stop("Parameter useGENO must be logical.", call. = FALSE)}
  if (! is.numeric(maf)) {stop("Parameter maf must be numeric.",call=FALSE)}  
  
  commandstr <- paste('plink','--bfile',  path2RefDir,
                              '--ld-snp-list', path2leadSNPList,
                              '--ld-window-kb', as.character(proxyWindow),
                              '--ld-window-r2',as.character(proxyCorr),
                              'r2')
  if (useMaf) {
    commandstr <-  paste(commandstr,'--maf',as.character(maf))
  }
  if (useHWE) {
    commandstr <-  paste(commandstr,'--hwe',as.character(hwe))
  }
  if (useGENO) {
    commandstr <-  paste(commandstr,'--geno',as.character(geno))
  }
  
  commandstr <-  paste(commandstr,'--out',as.character(path2Proxies))
  
  sys.run(commandstr)
  
}