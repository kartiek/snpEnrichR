#' Find proxies to SNPs using plink
#' 
#' ClumpSNPs is a wrapper to plink 1.9. The function computes proxy SNPs for the input snps 
#' and writes them to a file.
#' 
#' @param path2PlinkPrefix Plink bfile parameter. It is the path to reference directory including prefix of the plink reference files bed, bim and fam. For preprocessed files, see library(snpEnrichR).
#' @param path2leadSNPList Full path to the SNP list. 
#' @param ld_window_kb Plink parameter ld_window_kb denotes the maximum distance between LD buddies (default is 1000).
#' @param ld_window_r2 Plink parameter ld_window_r2 denoted the minimum correlation  between LD buddies (default is 0.8).
#' @param path2Proxies Path to directory where the result will be written. Must contain prefix for the output files.
#' @param ChainFilePath a UCSC chain format file to convert genome coordinates of the SNPs into genome build used by SNPsnap.  
#'    
#' @author Kari Nousiainen 
#' @export
#' 
#' @examples
#' findProxies(path2PlinkPrefix, snplist, path2Proxies)
#' @references Chang, C. C., Chow, C. C., Tellier, L. C., Vattikuti, S., Purcell, S. M., Lee, J. J. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. \emph{Gigascience}, 4(1), 7.
#' @references Purcell S.M., Chang C.C. \emph{PLINK 1.9} \url{ www.cog-genomics.org/plink/1.9/}.
#' @references Gaunt, T. R., Rodríguez, S., Day, I. N. (2007). Cubic exact solutions for the estimation of pairwise haplotype frequencies: implications for linkage disequilibrium analyses and a web tool'CubeX'. \emph{BMC bioinformatics}, 8(1), 428.
findProxies <- function(path2PlinkPrefix,snplist,ld_window_kb=1000,ld_window_r2=0.8,path2Proxies,ChainFilePath=NULL)
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

  if (! is.null(ChainFilePath)) {
    ch = import.chain(ChainFilePath)
    scs <- unlist(sapply(snplist,function(x) strsplit(x,':'),USE.NAMES = F))
    snpTable <- data.frame(cbind(scs[seq(1,length(scs),by=2)],
                                 scs[seq(2,length(scs),by=2)],
                                 scs[seq(2,length(scs),by=2)]
    ),stringsAsFactors = F)
    colnames(snpTable)<-c('chr','start','end')
    
    snpObjects <- sort(unique(makeGRangesFromDataFrame(snpTable, 
                                                       seqnames.field="chr",
                                                       keep.extra.columns=TRUE)))
    seqlevelsStyle(snpObjects) <- "UCSC"
    snpObjectsLiftedOver<-unlist(liftOver(snpObjects, ch))
    seqlevelsStyle(snpObjectsLiftedOver) <- "Ensembl"
    snpsNewCoords<-data.frame(snpObjectsLiftedOver)
    snplist <- apply(snpsNewCoords[,1:2],
                     1,
                     function(x) 
                       paste(trimws(x), collapse =":"))
  }
    
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