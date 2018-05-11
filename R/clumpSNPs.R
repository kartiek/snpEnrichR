#' Decipher a set of non-clumped SNPs from a list of SNPs
#' 
#' ClumpSNPs is a wrapper to plink 1.9. The function takes a character vector of SNPs and 
#' computes the clumped SNP in the vector. The output contains all SNPs without clumping and one 
#' representative for every clump. 
#'    
#' @param plinkPathPrefix Path to reference directory, must contain the prefix of the plink reference files bed, bim and fam. For preprocessed files, see library(snpEnrichR)
#' @param snplist Character vector of snps.
#' @param outputPath A character string is the Path to output file. 
#' @param clump_r2 R2 cutoff for clumping. 
#' @param clump_kb Distance cutoff for clumping.
#'
#' @return Returns a character vector of decorrelated SNPs.
#' @export 
#'
#' @examples
#' snpList <- clumpSNPs(plinkPathPrefix,snplist,outputPath,clump_r2,clump_kb)
#' @references Chang, C. C., Chow, C. C., Tellier, L. C., Vattikuti, S., Purcell, S. M., Lee, J. J. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. \emph{Gigascience}, 4(1), 7.
#' @references Purcell S.M., Chang C.C. \emph{PLINK 1.9} \url{ www.cog-genomics.org/plink/1.9/}.

clumpSNPs <- function(plinkPathPrefix,snplist,outputPath,clump_r2,clump_kb) {
    P=rep(0,length=length(snplist))
    snps <- cbind(SNP=snplist,P)
    
    print('SNPS read')
    tmpfilename <- tempfile(fileext = '.txt')

    write.table(snps,file=tmpfilename,quote=F,sep="\t",row.names = F)
    print('TMP file  created')
    cmdstr <- sprintf('plink --bfile %s --clump %s --clump-r2 %f --clump-kb %i --out %s',plinkPathPrefix,
                      tmpfilename,
                      clump_r2,
                      clump_kb,
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
  