#' Test for enrichment of SNPs within a given feature list.
#'
#' @param matchedSNPfile matched_snps.txt file from SNPsnap server
#' @param proxiesFile File containing proxies
#' @param featureObj Feature object for which encirhment has to be calculated
#'
#' @return Returns p.value
#' @export
#'
#' @examples
#' enrichTest()

enrichTest <- function(matchedSNPfile, proxiesFile, featureObj) {
  snpSets <- data.table::fread(matchedSNPfile, data.table = FALSE)
  allProx <- data.table::fread(proxiesFile,data.table = FALSE)
  getOvlp <- function(proxyList,randomSets,featureObj){
    proxyLocs <- proxyList %>% dplyr::filter(SNP_A %in% as.vector(unlist(randomSets))) %>% 
      dplyr::mutate(chr=paste0('chr',CHR_B)) %>% dplyr::select(chr,start=BP_B) %>% unique() %>% 
      with(GenomicRanges::GRanges(chr,IRanges(start = start,width = 1)))
    return(sum(IRanges::countOverlaps(proxyLocs,featureObj)))}
  ovList <- apply(snpSets,2,getOvlp,proxyList=allProx,featureObj=featureObj)
  disOv <- as.vector(unlist(ovList[1,1]))
  disP <- ovList[-1,] %>% dplyr::filter(ovList >= disOv) %>% nrow()
  return(disP/1000)
}