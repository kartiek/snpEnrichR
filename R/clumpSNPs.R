#' Title
#'
#' @param plinkPathPrefix 
#' @param snplist 
#' @param outputPath 
#' @param corrCoeff 
#' @param distance 
#'
#' @return Returns SNP list with clumped SNP buddies
#' @export
#'
#' @examples
clumpSNPs <- function(plinkPathPrefix,snplist,outputPath,corrCoeff,distance) {
    P=rep(0,length=length(snplist))
    snps <- cbind(SNP=snplist,P)
    
    print('SNPS read')
    tmpfilename <- tempfile(fileext = '.txt')

    write.table(snps,file=tmpfilename,quote=F,sep="\t",row.names = F)
    print('TMP file  created')
    cmdstr <- sprintf('plink --bfile %s --clump %s --clump-r2 %f --clump-kb %i --out %s',plinkPathPrefix,
                      tmpfilename,
                      corrCoeff,
                      distance,
                      tmpfilename)
    print(cmdstr)
    dd=system(cmdstr)

    res <- read.table(paste(tmpfilename,'clumped',sep="."),stringsAsFactors = F)
    colnames(res)<-res[1,]
    res <- res[-1,]
    nonclumpedSNPs <- res$SNP
    write.table(nonclumpedSNPs,file=outputPath,quote = F,row.names=F,col.names = F)
    return(nonclumpedSNPs)
    
}
  