#' Title
#'
#' @param plinkPathPrefix 
#' @param SNPsPath 
#' @param outputPath 
#' @param corrCoeff 
#' @param distance 
#'
#' @return
#' @export
#'
#' @examples
clumpSNPs <- function(plinkPathPrefix,SNPsPath,outputPath,corrCoeff,distance) {
  library(tidyverse)
  
    snps <- as.data.frame(read_tsv(SNPsPath,col_names="SNP",col_types = "c"))
    P=rep(0,length=nrow(snps))
    snps <- cbind(snps,P)
    
    print('SNPS read')
    tmpfilename <- tempfile(fileext = '.txt')

    write_tsv(snps,tmpfilename)
    print('TMP file  created')
    cmdstr <- sprintf('plink --bfile %s --clump %s --clump-r2 %f --clump-kb %i --out %s',plinkPathPrefix,
                      tmpfilename,
                      corrCoeff,
                      distance,
                      tmpfilename)
    print(cmdstr)
    system(cmdstr)
    
    res <- read.table(paste(tmpfilename,'clumped',sep="."),stringsAsFactors = F)
    colnames(res)<-res[1,]
    res <- res[-1,]
    nonclumpedSNPs <- as.data.frame(res$SNP)
    write_tsv(nonclumpedSNPs,outputPath,col_names = F)
    print(sprintf('Clumped SNPs written to %s',outputPath))  
    
}
  