#' analyzeEnrichment
#' 
#' Main function computes the overlaps between snps and genomic regions,
#' computes empirical p-values, performs Benjamini-Hochberg mutilple hypothesis
#' correction if needed, saves the results in a file. In addition, the function draws descriptive 
#' figures of the process.
#'
#' @param regionPath Path to file that contains the genomic regions. 
#' @param regionHeader Header of the genomics regions file.
#' @param SNPsnapPath Path to SNPsnap results.
#' @param numberOfRandomSNPsets The number of SNPs in a random SNP sets.
#' @param proxyPathPrefix  Path to the SNP file containing the file prefix (the extension is ld). The column SNP_A contains the names of the input SNPs. Columns CHR_B and BP_B contain the chromosome and the position of the SNP used in statistical hypothesis testing.   
#' @param traitShort Vector of abbreviated trait name.
#' @param traitsLong Vector of official name for traits.
#' @param genomicRegionsName The name describing the genomic regions.
#' @param cores The number of cores use in parallel computation. 
#' @param resDir Path to directory where the results are stored.
#' @param ChainFilePath a UCSC chain format file to map the genomic regions into hg19.
#' 
#' @return None
#' 
#' @export
#' 
#' @import GenomicFeatures 
#' @import dplyr
#' @import readr
#' @import parallel
#' 
#' @author Kari Nousiainen, Kartiek Kanduri
#' 
#' @examples
#'function(regionPath,SNPsnapPath,numberOfRandomSNPsets,proxyPathPrefix,traitShort,genomicRegionsName,cores,resDir)

analyzeEnrichment <-function(regionPath,
                             regionHeader=c('chr','start','end'),
                             SNPsnapPath,
                             numberOfRandomSNPsets,
                             proxyPathPrefix,
                             traitShort,
                             genomicRegionsName,
                             cores,
                             resDir,
                             traitsLong=NULL,
                             ChainFilePath=NULL)
{
  library(GenomicFeatures)
  library(readr)
  library(tidyr)
  library(dplyr)
  library(parallel)
  library(ggplot2)

  
  regions <- as.data.frame(read_tsv(regionPath))
  colnames(regions )[1:3] <-regionHeader

  f2 <- makeGRangesFromDataFrame(regions, 
                                 keep.extra.columns=FALSE,
                                 start.field="start",
                                 end.field="end",
                                 starts.in.df.are.0based=FALSE)
  if (! is.null(ChainFilePath)) {
    seqlevelsStyle(f2) <- "UCSC"
    ch = import.chain(ChainFilePath)
    f2<-unlist(liftOver(f2, ch))
  }
  f2=reduce(f2)
  
  # Read SNP data
  getSNPDat <- function(x){
  nc=max(count.fields(file.path(paste(SNPsnapPath,x,sep=''),"matched_snps.txt"),sep="\t"))
  nc=min(numberOfRandomSNPsets+1,nc)
  snpDat <- read.table(file.path(paste(SNPsnapPath,x,sep=''),'matched_snps.txt'),header=T,stringsAsFactors=F) %>% dplyr::select(1:nc)

  return(snpDat)}

  # Read SNP proxies
  readProx <- function(x){
    return(read.table(paste0(proxyPathPrefix,x,'.ld'),sep = '',header = T))
  }

# Get proxies for each set and their respective overlaps
  getOl <- function(w,y,z){
    y1 <- w %>% dplyr::filter(SNP_A %in% as.vector(unlist(y))) %>% dplyr::select(CHR_B,BP_B) %>% unique()
    colnames(y1) <- c('chr','start')
    y1$chr <- paste0('chr',y1$chr)
    f1 <- with(y1,GRanges(chr,IRanges(start = start,width = 1)))
    return(sum(countOverlaps(f1,z)))}

# Analyse the data and write results
  analyzeOVs <- function(x){
    snpSets <- getSNPDat(x)
    allProx <- readProx(x)
    ovList <- apply(snpSets,2,getOl,w=allProx,z=f2) 
    write_tsv(as.data.frame(ovList),file.path(resDir,paste0('Overlaps_',x,'.txt')))}


  disRes <-  mclapply(traitShort,analyzeOVs,mc.cores = cores)

  getPVal <- function(x){
    disFile <- read.table(file.path(resDir,paste0('Overlaps_',x,'.txt')),header=T)
    disOv <- as.vector(unlist(disFile[1,1]))
    disP <- disFile[-1,,drop=F] %>% 
    filter(ovList >= disOv) %>% 
    nrow()
    return(disP/1000)}

  disPval <- data_frame(disease=traitShort,pval=p.adjust(sapply(traitShort,getPVal),method = 'BH'))
  ef <- 'disease_enrichment.txt'
  if(!file.exists(file.path(resDir,ef))) {
    write_tsv(disPval,file.path(resDir,ef))
  } else {
    write_tsv(disPval,file.path(resDir,ef),append=T,col_names=F)
  }
  # Plot the background dists
  allBkg <- bind_cols(lapply(file.path(resDir,paste0('Overlaps_',traitShort,'.txt')),read_tsv))
  colnames(allBkg) <- traitShort
  allBD <- allBkg[1,]
  disVal <- gather(allBD,trait,val)
  allBkg <- allBkg[-1,]
  bkg <- gather(allBkg,dis,val)

  pvb <- read_tsv(file.path(resDir,ef))
  colnames(pvb) <- c('dis','pval')
  if (nrow(pvb)>1) {
    png(file.path(resDir,'background_distributions.png'),width = 16,height = 9,units = 'in',res = 300)
  } else {
    png(file.path(resDir,paste0(traitShort,'_','background_distributions.png')),width = 16,height = 9,units = 'in',res = 300)
  } 
  
  print(ggplot(bkg,aes(val)) + geom_histogram(binwidth=1) +
          facet_wrap(~dis,ncol = 4,scales = 'free') + theme_bw(base_size = 16) +
          geom_vline(data=disVal,aes(xintercept = val),linetype='dotted') +
          geom_text(data = pvb,aes(label=paste("adj.p.value =", round(pval,3))),inherit.aes = T,x=Inf,y=Inf,vjust=1.5,hjust=1.2) +
          labs(title='Background distributions of diseases',subtitle='Dotted line indicates the actual value',x=NULL,y=NULL))
  dev.off()

  # Plot the number of SNPs for diseases
  if(is.null(traitsLong)) { 
      traitsLong <- traitShort
  }
  disVal <- disVal[-(length(traitsLong)+1),]
  disVal$dis <-traitsLong
  png(file.path(resDir,'autImmDisSNPs.png'),width=9,height=9,res=300,units = 'in')
      ggplot(disVal,aes(dis,val)) + geom_bar(stat='identity') + 
      ggthemes::theme_tufte(base_size=18) + coord_flip() + 
      labs(title=sprintf('Lead SNPs and their proxies from \n %s overlapping %s genomic regions', traitShort,genomicRegionsName),y='Number of SNPs',x=NULL) + 
      theme(axis.text = element_text(face = 'bold',size = 18))
  dev.off()



}