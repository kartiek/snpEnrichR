library(tidyverse)
removeHighlyCorrelatedSNPs <- function(proxyPathPrefix,SNPsPath,outputPath,corrCoeff,distance) {
    snps <- as.data.frame(read_tsv(SNPsPath,col_names="SNP",col_types = "c"))
    P=rep(0,length=nrow(snps))
    snps <- cbind(snps,P)
    tmpfilename <- paste(SNPsPath,"tmp.assoc",sep=".")
    tmpfilename2 <- paste(tmpfilename,"clumped",sep=".")
    write_tsv(snps,tmpfilename)
    cmdstr <- sprintf('plink --bfile %s --clump %s --clump-r2 %f --clump-kb %i --out %s',proxyPathPrefix,
                      tmpfilename,
                      corrCoeff,
                      distance,
                      tmpfilename2)
    system(cmdstr)
    res <- read.table(tmpfilename2,stringsAsFactors = F)
    colnames(res)<-res[1,]
    res <- res[-1,]
    nonclumpedSNPs <- as.data.frame(res$SNP)
    write_tsv(nonclumpedSNPs,outputPath,col_names = F)
  
    
}
  