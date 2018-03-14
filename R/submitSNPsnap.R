#' Submit list of SNPs to SNPsnap server
#'
#' @param snplist List of SNPs
#' @param super_population Super population to use; defaults to 'EUR'
#' @param distance_type Distance type; default is 'ld'
#' @param distance Distance to use
#' @param max_freq_deviation Percentage deviation from query SNP minor allele frequency
#' @param max_genes_count_deviation Gene density deviation from query SNP
#' @param max_distance_deviation Gene distance deviation from query SNP
#' @param max_ld_buddy_count_deviation Deviation of number of LD buddies near query SNP 
#' @param ld_buddy_cutoff Cutoff of R2 while calculation LD buddies
#' @param N_sample_sets Number of SNP sets requested
#' @param annotate_matched a logical value to annotate matched SNPs
#' @param annotate_input a logical value to annotate input SNPs
#' @param clump_input a logical value to report clumping of input SNPs
#' @param clump_r2 R2 cutoff for clumping
#' @param clump_kb Distance cutoff for clumping
#' @param exclude_input_SNPs a logical value to exclude input SNPs in matched SNPs
#' @param exclude_HLA_SNPs a logical value to exclude HLA SNPs
#' @param email_address Your email address
#' @param job_name Job name
#'
#' @return Returns a URL from which results can be downloaded
#' @export
#' @import RSelenium
#' @importFrom readr write_tsv
#' @examples
#' submitSNPsnap(snplist=snps,email_address='kartiek.kanduri@aalto.fi')

submitSNPsnap <- function(snplist, super_population = c('EUR','EAS','WAFR'),
                          distance_type = c('ld','kb'), distance, max_freq_deviation = seq(1,50,1),
                          max_genes_count_deviation, max_distance_deviation,
                          max_ld_buddy_count_deviation, ld_buddy_cutoff = seq(0.1,0.9,0.1),
                          N_sample_sets, annotate_matched = FALSE, annotate_input = FALSE,
                          clump_input = TRUE, clump_r2 = seq(0.1,0.9,0.1), clump_kb = seq(100,1000,100),
                          exclude_input_SNPs = TRUE, exclude_HLA_SNPs = TRUE,
                          email_address, job_name){
  
  rD <- rsDriver(verbose = FALSE, browser = 'phantomjs')
  remDr <- rD$client
  
  remDr$navigate(snpSnapURL)
  
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
