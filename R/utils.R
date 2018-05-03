
# EBI GWAS URL
gwasURL <- 'https://www.ebi.ac.uk/gwas/api/search/downloads?q=&pvalfilter=&orfilter=&betafilter=&datefilter=&genomicfilter=&traitfilter%5B%5D=&dateaddedfilter=&facet=association'

#SNPsnap URL
snpSnapURL <- 'https://data.broadinstitute.org/mpg/snpsnap/match_snps.html'

# This following split_to_df function is originally from rsnps package but modified to my needs.
split_to_df <- function(x, sep, fixed=FALSE, perl=TRUE, useBytes=FALSE, names=NULL) {
  
  x <- as.character(x)
  
  if( fixed ) {
    perl <- FALSE
  }
  
  tmp <- strsplit( x, sep, fixed=fixed, perl=perl, useBytes=useBytes )
  if( length( unique( unlist( lapply( tmp, length ) ) ) ) > 1 ) {
    stop("non-equal lengths for each entry of x post-splitting")
  }
  tmp <- unlist( tmp )
  tmp <- as.data.frame( 
    matrix( tmp, ncol = (length(tmp) / length(x)), byrow=TRUE ),
    stringsAsFactors=FALSE, optional=TRUE 
  )
  
  if( !is.null(names) ) {
    names(tmp) <- names
  } else {
    names(tmp) <- tmp[1,]
    tmp <- tmp[-1,]
  }
  
  return(tmp)
}

# Build query for SNPsnap submission
buildQuery <- function(snplist, super_population = c('EUR','EAS','WAFR'),
                       distance_type = c('ld','kb'), distance, max_freq_deviation = seq(1,50,1),
                       max_genes_count_deviation, max_distance_deviation,
                       max_ld_buddy_count_deviation, ld_buddy_cutoff = seq(0.1,0.9,0.1),
                       N_sample_sets, annotate_matched = FALSE, annotate_input = FALSE,
                       clump_input = TRUE, clump_r2 = seq(0.1,0.9,0.1), clump_kb = seq(100,1000,100),
                       exclude_input_SNPs = TRUE, exclude_HLA_SNPs = TRUE,
                       email_address, job_name){
  
  snplist_fileupload <- tempfile(fileext = '.txt')
  write_tsv(as.data.frame(snplist),snplist_fileupload)
  
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
 
  argsList <- list(snplist_fileupload = snplist_fileupload,
                   super_population = super_population,
                   distance_type = distance_type,
                   distance = distance,
                   max_freq_deviation = max_freq_deviation,
                   max_genes_count_deviation = max_genes_count_deviation,
                   max_distance_deviation = max_distance_deviation,
                   max_ld_buddy_count_deviation = max_ld_buddy_count_deviation,
                   ld_buddy_cutoff = ld_buddy_cutoff,
                   N_sample_sets = N_sample_sets,
                   annotate_matched = annotate_matched,
                   annotate_input = annotate_input,
                   clump_input = clump_input,
                   clump_r2 = clump_r2,
                   clump_kb = clump_kb,
                   exclude_input_SNPs = exclude_input_SNPs,
                   exclude_HLA_SNPs = exclude_HLA_SNPs,
                   email_address = email_address,
                   job_name = job_name)
  
  return(argsList)

}

# Function to fill data into SNPsnap submission form
sendData <- function(driver, argName) {
  driver$findElement(using = 'name', value = argName)$clearElement()
  driver$findElement(using = 'name', value = argName)$sendKeysToElement(list(as.character(argsList$argName)))
}