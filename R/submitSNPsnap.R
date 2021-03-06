#' Submit list of SNPs to SNPsnap server
#' 
#' SubmitSNPsnap submits a query to SNPsnap server for creating random sets of SNPs matching input SNPs and 
#' the criteria given by the user.
#'
#' @param snplist List of SNPs.
#' @param super_population Super population to use; defaults to 'EUR'.
#' @param distance_type Distance type; default is 'ld'.
#' @param distance Distance to use.
#' @param max_freq_deviation Percentage deviation from query SNP minor allele frequency.
#' @param max_genes_count_deviation Gene density deviation from query SNP.
#' @param max_distance_deviation Gene distance deviation from query SNP.
#' @param max_ld_buddy_count_deviation Deviation of number of LD buddies near query SNP. 
#' @param ld_buddy_cutoff Cutoff of R2 while calculation LD buddies.
#' @param N_sample_sets Number of SNP sets requested.
#' @param annotate_matched a logical value to annotate matched SNPs.
#' @param annotate_input a logical value to annotate input SNPs.
#' @param clump_input a logical value to report clumping of input SNPs.
#' @param clump_r2 R2 cutoff for clumping.
#' @param clump_kb Distance cutoff for clumping.
#' @param exclude_input_SNPs a logical value to exclude input SNPs in matched SNPs.
#' @param exclude_HLA_SNPs a logical value to exclude HLA SNPs.
#' @param email_address Your email address.
#' @param ChainFilePath a UCSC chain format file to convert genome coordinates of the SNPs into genome build used by SNPsnap.  
#' @param job_name Job name. SNPsnap results are written in a zip file called SNPsnap_jobname.
#'
#' @return Returns a URL from which results can be downloaded
#' @author Kartiek Kanduri, Kari Nousiainen
#' @export
#' @import RSelenium
#' @import rtracklayer
#' @importFrom readr write_tsv
#' @examples
#' submitSNPsnap(snplist=snps,email_address='kartiek.kanduri@aalto.fi')
#' @references Pers, T. H., Timshel, P., Hirschhorn, J. N. (2014). SNPsnap: a Web-based tool for identification and annotation of matched SNPs. \emph{Bioinformatics}, 31(3), 418-420.
submitSNPsnap <- function(snplist, super_population = c('EUR','EAS','WAFR'),
                          distance_type = c('ld','kb'), distance, max_freq_deviation = seq(1,50,1),
                          max_genes_count_deviation, max_distance_deviation,
                          max_ld_buddy_count_deviation, ld_buddy_cutoff = seq(0.1,0.9,0.1),
                          N_sample_sets, annotate_matched = FALSE, annotate_input = FALSE,
                          clump_input = TRUE, clump_r2 = seq(0.1,0.9,0.1), clump_kb = seq(100,1000,100),
                          exclude_input_SNPs = TRUE, exclude_HLA_SNPs = TRUE,
                          email_address, job_name,ChainFilePath = NULL) {

  library(RSelenium)
  library(rtracklayer)
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
  tFile <- tempfile(fileext = '.txt')
  readr::write_tsv(as.data.frame(snplist),tFile)
  if(missing(super_population)){super_population <- 'EUR'}
  else{ super_population <- match.arg(super_population) }
  if (missing(distance_type)){distance_type <- 'ld'}
  else{ distance_type <- match.arg(distance_type) }
  valid_distances <- distance_type
  valid_distances <- switch (valid_distances, ld = seq(0.1,0.9,0.1), kb = seq(100,1000,100))
  if ( !(distance %in% valid_distances) ) {
    stop("invalid distance. distance must be one of: ",
         paste( valid_distances, collapse = ", "), call. = FALSE)
  }
  if( missing(max_freq_deviation) ){ max_freq_deviation <- 5 }
  else { max_freq_deviation <- match.arg(max_freq_deviation)}
  if(missing(max_genes_count_deviation)){max_genes_count_deviation <- 50}
  else if ( max_genes_count_deviation < 1 ) {
    stop("invalid max_genes_count_deviation. max_genes_count_deviation must be greater than or equal to 1", call. = FALSE)
  }
  if(missing(max_distance_deviation)){max_distance_deviation <- 50}
  else if ( max_distance_deviation < 1 ) {
    stop("invalid max_distance_deviation. max_distance_deviation must be greater than or equal to 1", call. = FALSE)
  }
  if(missing(max_ld_buddy_count_deviation)){max_ld_buddy_count_deviation <- 50}
  else if ( max_ld_buddy_count_deviation < 1 ) {
    stop("invalid max_ld_buddy_count_deviation. max_ld_buddy_count_deviation must be greater than or equal to 1", call. = FALSE)
  }
  valid_lds <- seq(0.1,0.9,0.1)
  if(missing(ld_buddy_cutoff)){ld_buddy_cutoff <- 0.5}
  else if ( !(ld_buddy_cutoff %in%  valid_lds) ) {
    stop("invalid ld. Value must be one of: ",
         paste( valid_lds, collapse = ", "), call. = FALSE)
  }
  else{ld_buddy_cutoff <- match.arg(ld_buddy_cutoff)}
  if(missing(N_sample_sets)){N_sample_sets <- 10000}
  else if ( N_sample_sets < 1 || N_sample_sets > 20000) {
    stop("invalid N_sample_sets. N_sample_sets must be between 1 and 20000", call. = FALSE)
  }
  if(missing(annotate_matched)){annotate_matched <- FALSE}
  if(missing(annotate_input)){annotate_input <- FALSE}
  if(missing(clump_input)){clump_input <- TRUE}
  if ( clump_input == TRUE ){
    if ( !(clump_r2 %in%  valid_lds) ) {
      stop("invalid clumping ld. Value must be one of: ",
           paste( valid_lds, collapse = ", "), call. = FALSE)
    }
    valid_kbs = seq(100,1000,100)
    if ( !(clump_kb %in%  valid_kbs) ) {
      stop("invalid clumping distance. Value must be one of: ",
           paste( valid_kbs, collapse = ", "), call. = FALSE)
    }
    clump_r2 <- match.arg(clump_r2)
    clump_kb <- match.arg(clump_kb)
  }
  if(missing(exclude_input_SNPs)){exclude_input_SNPs <- TRUE}
  if(missing(exclude_HLA_SNPs)){exclude_HLA_SNPs <- TRUE}
  if(missing(job_name)){job_name <- 'testJob'}
  if(! is.character(job_name))  {stop("Parameter job_name should be a string.", call. = FALSE)}
  if(missing(email_address)){
    stop("Please provide an email address", call. = FALSE) }
  else if(! is.character(email_address))  {stop("Parameter email_address should be a string.", call. = FALSE)}
  rD <- rsDriver(verbose = FALSE)#, browser = 'phantomjs')
  #rD <- rsDriver(verbose = FALSE, browser = 'phantomjs')
  remDr <- rD$client
  
  remDr$navigate(snpSnapURL)
  snplist_fileupload<-tFile
  aplyArgs <- c(snplist_fileupload, super_population, max_freq_deviation, max_genes_count_deviation,
                max_distance_deviation, max_ld_buddy_count_deviation, ld_buddy_cutoff, 
                N_sample_sets, job_name, email_address)
  
  
  
  remDr$findElement(using = 'name', value = "snplist_fileupload")$sendKeysToElement(list(tFile))
  remDr$findElement(using = 'name', value = "super_population")$sendKeysToElement(list(super_population))
  if ( distance_type == 'kb' ){
    remDr$findElement(using = 'id', value = 'kb_distance_type')$clickElement()
    remDr$findElement(using = 'id', value = "kb_distance_cutoff")$sendKeysToElement(list(as.character(distance)))
  }
  else if ( distance_type == 'ld' ){
    remDr$findElement(using = 'id', value = 'ld_distance_type')$clickElement()
    remDr$findElement(using = 'id', value = "ld_distance_cutoff")$sendKeysToElement(list(as.character(distance)))
  }
  remDr$findElement(using = 'name', value = "max_freq_deviation")$clearElement()
  remDr$findElement(using = 'name', value = "max_freq_deviation")$sendKeysToElement(list(as.character(max_freq_deviation)))
  remDr$findElement(using = 'name', value = "max_genes_count_deviation")$clearElement()
  remDr$findElement(using = 'name', value = "max_genes_count_deviation")$sendKeysToElement(list(as.character(max_genes_count_deviation)))
  remDr$findElement(using = 'name', value = "max_distance_deviation")$clearElement()
  remDr$findElement(using = 'name', value = "max_distance_deviation")$sendKeysToElement(list(as.character(max_distance_deviation)))
  remDr$findElement(using = 'name', value = "max_ld_buddy_count_deviation")$clearElement()
  remDr$findElement(using = 'name', value = "max_ld_buddy_count_deviation")$sendKeysToElement(list(as.character(max_ld_buddy_count_deviation)))
  
  remDr$findElement(using = 'name', value = "ld_buddy_cutoff")$sendKeysToElement(list(as.character(ld_buddy_cutoff)))
  remDr$findElement(using = 'name', value = "N_sample_sets")$clearElement()
  remDr$findElement(using = 'name', value = "N_sample_sets")$sendKeysToElement(list(as.character(N_sample_sets)))
  if(annotate_matched == TRUE){
    if(remDr$findElement(using = 'id', value = 'set_file')$isElementSelected() == FALSE){
      remDr$findElement(using = 'id', value = 'set_file')$clickElement()
    }}
  if(annotate_input == TRUE){
    if(remDr$findElement(using = 'id', value = 'annotate')$isElementSelected() == FALSE){
      remDr$findElement(using = 'id', value = 'annotate')$clickElement()
    }}
  if(clump_input == TRUE){
    if(remDr$findElement(using = 'id', value = 'clump')$isElementSelected() == FALSE){
      remDr$findElement(using = 'id', value = 'clump')$clickElement()
      remDr$findElement(using = 'name', value = "clump_r2")$sendKeysToElement(list(as.character(clump_r2)))
      remDr$findElement(using = 'name', value = "clump_kb")$clearElement()
      remDr$findElement(using = 'name', value = "clump_kb")$sendKeysToElement(list(as.character(clump_kb)))
    }}
  if(exclude_input_SNPs == TRUE){
    if(remDr$findElement(using = 'id', value = 'exclude_input_SNPs')$isElementSelected() == FALSE){
      remDr$findElement(using = 'id', value = 'exclude_input_SNPs')$clickElement()
    }}
  if(exclude_HLA_SNPs == TRUE){
    if(remDr$findElement(using = 'id', value = 'exclude_HLA_SNPs')$isElementSelected() == FALSE){
    remDr$findElement(using = 'id', value = 'exclude_HLA_SNPs')$clickElement()
    }}
  remDr$findElement(using = 'name', value = "job_name")$clearElement()
  remDr$findElement(using = 'name', value = "job_name")$sendKeysToElement(list(job_name))
  remDr$findElement(using = 'name', value = "email_address")$sendKeysToElement(list(email_address))
  remDr$findElement(using = 'class', value = 'btn-success')$clickElement()
  webElem <- remDr$findElement(using = 'class', value = 'btn-success')
  
  urlOut <- webElem$getElementAttribute('href')[[1]]
  
  message('Results can be downloaded from ',webElem$getElementAttribute('href')[[1]])
  
  remDr$close()
  rD$server$stop()
  
  return(urlOut)
}
