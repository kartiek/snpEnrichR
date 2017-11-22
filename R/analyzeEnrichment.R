
#' Title
#'
#' @param regionPath 
#' @param regionHeader 
#' @param SNPsnapPath 
#' @param numberOfRandomSNPsets 
#' @param proxyPathPrefix 
#' @param traitShort 
#' @param traitsLong 
#' @param genomicRegionsName 
#' @param cores
#' @param resDir 
#'
#' @return
#' @export
#'
#' @examples
analyzeEnrichment <-function(regionPath,regionHeader=c('chr','start','end'),SNPsnapPath,numberOfRandomSNPsets,LSProxyPathPrefix,BGProxyPathPrefix,traitShort,genomicRegionsName,cores,resDir,traitsLong=NULL)
{
  library(GenomicFeatures)
  library(tidyverse)
  library(parallel)

  
  regions <- as.data.frame(read.table(regionPath))
  colnames(regions ) <-regionHeader
  f2 <- with(regions ,GRanges(chr,IRanges(start = start,end = end)))

  # Read SNP data
  getSNPDat <- function(x){
  nc=max(count.fields(file.path(paste(SNPsnapPath,x,sep=''),"matched_snps.txt"),sep="\t"))
  nc=min(numberOfRandomSNPsets+1,nc)
  snpDat <- read.table(file.path(paste(SNPsnapPath,x,sep=''),'matched_snps.txt'),header=T,stringsAsFactors=F) %>% dplyr::select(1:nc)

  return(snpDat)}

  # Read SNP proxies
  readProx <- function(x,proxyPathPrefix){
    return(read.table(file.path(proxyPathPrefix,paste0(x,'.ld')),sep = '',header = T,stringsAsFactors=F))
    }

# Get proxies for each set and their respective overlaps
  getOl <- function(w,y,z){
    y1 <- w %>% dplyr::filter(SNP_A %in% as.vector(unlist(y))) %>% dplyr::select(4,5) %>% unique()
    colnames(y1) <- c('chr','start')
    y1$chr <- paste0('chr',y1$chr)
    f1 <- with(y1,GRanges(chr,IRanges(start = start,width = 1)))
    return(sum(countOverlaps(f1,z)))}

# Analyse the data and write results
  analyzeOVs <- function(x){
    snpSets <- getSNPDat(x)
    LeadSNPProx <- readProx(x,LSProxyPathPrefix)
    BGPSNPProx <- readProx(x,BGProxyPathPrefix)
    ovList1 <- apply(snpSets[,1,drop=F],2,getOl,w=LeadSNPProx,z=f2)
    ovList2 <- apply(snpSets[,2:(numberOfRandomSNPsets+1),drop=F],2,getOl,w=BGPSNPProx,z=f2)
    ovList=data.frame(ovList=c(ovList1,ovList2))
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
  write_tsv(disPval,file.path(resDir,'disease_enrichment.txt'))

  # Plot the background dists
  allBkg <- bind_cols(lapply(file.path(resDir,paste0('Overlaps_',traitShort,'.txt')),read_tsv))
  colnames(allBkg) <- traitShort
  allBD <- allBkg[1,]
  disVal <- gather(allBD,trait,val)
  allBkg <- allBkg[-1,]
  bkg <- gather(allBkg,dis,val)

  pvb <- read_tsv(file.path(resDir,'disease_enrichment.txt'))
  colnames(pvb) <- c('dis','pval')
  png(file.path(resDir,'background_distributions.png'),width = 16,height = 9,units = 'in',res = 300)
  print(ggplot(bkg,aes(val)) + geom_histogram(binwidth=1) +
          facet_wrap(~dis,ncol = 4,scales = 'free') + theme_bw(base_size = 16) +
          geom_vline(data=disVal,aes(xintercept = val),linetype='dotted') +
          geom_text(data = pvb,aes(label=paste("adj.p.value =", round(pval,3))),inherit.aes = T,x=Inf,y=Inf,vjust=1.5,hjust=1.2) +
          labs(title='Background distributions of diseases',subtitle='Dotted line indicates the actual value',x=NULL,y=NULL))
  dev.off()

  # Plot the number of SNPs for diseases
  if(is.null(traitsLong)) { 
      traitsLong <- traitShort#c('Ankylosing spondylitis','Celiac disease','Crohn\'s disease','Immunoglobulin A','Multiple sclerosis','Primary biliary cirrhosis','Psoriasis','Rheumatoid arthritis','Systemic lupus erythematosus','Type 1 diabetes','Ulcerative colitis')
  }
  disVal <- disVal[-(length(traitsLong)+1),]
  disVal$dis <-traitsLong
  png(file.path(resDir,'autImmDisSNPs.png'),width=9,height=9,res=300,units = 'in')
      ggplot(disVal,aes(dis,val)) + geom_bar(stat='identity') + 
      ggthemes::theme_tufte(base_size=18) + coord_flip() + 
      labs(title=sprintf('Lead SNPs and their proxies from \n %s overlapping %s genomicregions', traitShort,genomicRegionsName),y='Number of SNPs',x=NULL) + 
      theme(axis.text = element_text(face = 'bold',size = 18))
  dev.off()



}