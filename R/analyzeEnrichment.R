library(GenomicFeatures)
library(tidyverse)
library(parallel)

analyzeEnrichment <-function(PeakPath,peakHeader,SNPsnapPath,proxyPathPrefix,traitShort,traitsLong,traitName,genomicRegionsName,cores)
{                             
peaks <- as.data.frame(read_tsv(PeakPath))
colnames(peaks ) <-peaskHeader
f2 <- with(peaks ,GRanges(chr,IRanges(start = start,end = end),name=name))

# Read SNP data
getSNPDat <- function(x){
  snpDat <- read_tsv(file.path(paste(SNPsnapPath,x,sep=''),'matched_snps.txt'),col_types = paste0(rep('c',1001),collapse = '')) %>% 
    dplyr::select(1:1001)
  return(snpDat)}

# Read SNP proxies
readProx <- function(x){
  return(read.table(paste0(proxyPathPrefix,x,'.ld'),sep = '',header = T))
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
  allProx <- readProx(x)
  ovList <- apply(snpSets,2,getOl,w=allProx,z=f2)
  write_tsv(as.data.frame(ovList),paste0('Overlaps_',x,'.txt'))}


disRes <-  mclapply(traits,analyzeOVs,mc.cores = cores)
stopCluster(cluster)

getPVal <- function(x){
  disFile <- read_tsv(paste0('Overlaps_',x,'.txt'))
  disOv <- as.vector(unlist(disFile[1,1]))
  disP <- disFile[-1,] %>% 
    filter(ovList >= disOv) %>% 
    nrow()
  return(disP/1000)}
disPval <- data_frame(disease=traits,pval=p.adjust(sapply(traits,getPVal),method = 'BH'))
write_tsv(disPval,'disease_enrichment.txt')

# Plot the background dists
allBkg <- bind_cols(lapply(paste0('Overlaps_',traits,'.txt'),read_tsv))
colnames(allBkg) <- traits
allBD <- allBkg[1,]
disVal <- gather(allBD,trait,val)
allBkg <- allBkg[-1,]
bkg <- gather(allBkg,dis,val)

pvb <- read_tsv('disease_enrichment.txt')
colnames(pvb) <- c('dis','pval')
png('background_distributions.png',width = 16,height = 9,units = 'in',res = 300)
print(ggplot(bkg,aes(val)) + geom_density() +
        facet_wrap(~dis,ncol = 4,scales = 'free') + theme_bw(base_size = 16) +
        geom_vline(data=disVal,aes(xintercept = val),linetype='dotted') +
        geom_text(data = pvb,aes(label=paste("adj.p.value =", round(pval,3))),inherit.aes = T,x=Inf,y=Inf,vjust=1.5,hjust=1.2) +
        labs(title='Background distributions of diseases',subtitle='Dotted line indicates the actual value',x=NULL,y=NULL))
dev.off()

# Plot the number of SNPs for diseases
traitsLong <- c('Ankylosing spondylitis','Celiac disease','Crohn\'s disease','Immunoglobulin A','Multiple sclerosis','Primary biliary cirrhosis','Psoriasis','Rheumatoid arthritis','Systemic lupus erythematosus','Type 1 diabetes','Ulcerative colitis')
disVal <- disVal[-(length(traitsLong)+1),]
disVal$dis <-traitsLong
png('autImmDisSNPs.png',width=9,height=9,res=300,units = 'in')
ggplot(disVal,aes(dis,val)) + geom_bar(stat='identity') + 
  ggthemes::theme_tufte(base_size=18) + coord_flip() + 
  labs(title=sprintf('Lead SNPs and their proxies from \n %s overlapping %s peaks', traitName,genomicRegionsName),y='Number of SNPs',x=NULL) + 
  theme(axis.text = element_text(face = 'bold',size = 18))
dev.off()



}