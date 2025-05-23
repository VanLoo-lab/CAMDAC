##  Get differential methylation calls
##  Latest update: 04/05/2021
##  Version 1.0.0
##  Author: Elizabeth Larose Cadieux

#' Perform differential methylation analysis on deconvolute tumour methylation rates
#' 
#' \code{get_differential_methylation}
#'
#' @param patient_id Character variable containting the patient id number
#' @param sample_id Character variable with the tumour sample_id
#' @param sex Character variable with  the patient expressed as "XX" for
#' female or "XY" for male.
#' @param normal_origin_proxy_id Character variable with the sample ID
#' of the normal to be used as a proxy for the tumour cell of origin in
#' @param path Character path variable pointing to the desired working
#' directory. This is where the output will be stored.
#' @param path_to_CAMDAC Character variable containting the path to the CAMDAC
#' directory including dir name (e.g. "/path/to/CAMDAC/").
#' @param build Character variable corresponding to the reference genome
#' used for alignment. CAMDAC is compatible with "hg19", "hg38", "GRCH37","GRCH38".
#' @param effect_size Numerical containting the minimum tumour-normal methylation
#' difference (default is 0.2)
#' @param prob Numerical  value representing the threshold for statistically
#' significant DMP (default is p=0.99)
#' @param min_DMP_counts_in_DMR Numerical value representing the number of
#' DMPs required in a DMR
#' @param min_consec_DMP_in_DMR Numerical value representing the number of
#' consecutive DMPs required in a DMR
#' @param n_cores Numerical value correspdonding to the number of cores
#' for parallel processing
#' @param reseg Logical value should be set to FALSE. Multi-sample re-segmentation of
#' the copy number profiles will be available in future versions of CAMDAC.
#' @param bulk Default is FALSE unless you want bulk DMP/DMR calls in addition
#' to CAMDAC pure tumour differential methylation analysis
#'
#' Note:
#' #' Annotation include:
#' CGI (including shore and shelves)
#' gene body (intragenic, 5UTR, 3UTR, intron, exon)
#' promoter (2kb upstream and 500 downstream any UCSC annotated gene)
#' enhancer (vista and FANTOM5 annotation)
#'  
#' @return Biologically significant DMPs, DMRs

get_differential_methylation <-
    function(patient_id,sample_id,sex,normal_origin_proxy_id,
             path,path_to_CAMDAC,build,
             effect_size=0.2,prob=0.99,
             min_DMP_counts_in_DMR=5,min_consec_DMP_in_DMR=4,
             n_cores, reseg=FALSE, bulk=FALSE){
  
  control_id=normal_origin_proxy_id
  if(sample_id==control_id){
      stop(paste("Differential methylation will be obtained between tumour sample X",
                 "and your proxy for the normal cell of origin.",
                 "As such, you cannot set sample_id to your normal cell of origin sample ID.",
                 sep="\n"))
  }
  if(parallel::detectCores()<n_cores){
      warning(paste0(n_cores, " cores selected but only ",
      parallel::detectCores(), " detected on machine."))
      
  }
  
  # set build to to assembly version disregarging UCSC/Ensembl
  if(build=="GRCH37"){build="hg19"}
  if(build=="GRCH38"){build="hg38"}
  
  # set output path
  path_patient <- file.path(path, patient_id)
  path_output <-paste0(path_patient, "/Methylation/",sample_id,"/")
  if(reseg == TRUE){ path_output <- paste0(path_output, "/re_segmented/") }
    
  # load purified tumour data.table
  load(paste0(path_output, "purified_tumour", ".RData"))
  dt <- dt_purified_tumour
  
  # Ensure dt is a data.table object
  if(class(dt)[1]!="data.table"){
  dt <- data.table(dt)
  }

  # remove columns with normal infiltrates info and keep only normal contaminants
  if(sum(grepl("^M_n_o$", colnames(dt)))==1){
      dt <- dt[, .SD, .SDcols=c(1:ncol(dt))[!grepl("_n_i", colnames(dt))]]
      colnames(dt)[grepl("_n_o", colnames(dt))] <-
                gsub("_o", "", colnames(dt)[grepl("_n_o", colnames(dt))])
  }
 
  # annotate tumour-normal methylation difference
  if(bulk==TRUE){ dt[, m_diff_bn := (m_b-m_n)] }# bulk tumour - normal
  dt[, m_diff_tn := (m_t_corr-m_n)] # CAMDAC pure tumour - normal

  # get DMPs
  dt <- get_DMPs(path=path_output, # results path to save DMP stats
                 patient_id=patient_id, sample_id=sample_id, # for DMP stats filename
                 df=dt, # methylation info for DMP calling
                 n_cores=n_cores) # number of cores for parallel processing
  
  # identify biologically significant DMPs in the bulk
  if(bulk==TRUE){ 
  dt[, DMP_b := ifelse(prob_DMP >= prob & m_diff_bn >= effect_size, "hyper",
                      ifelse(prob_DMP >= prob & m_diff_bn <= (-effect_size), "hypo", NA))]
  }
    
  # identify biologically significant DMPs in the pure tumour
  dt[, DMP_t := ifelse(prob_DMP >= prob & m_diff_tn >= effect_size, "hyper",
                      ifelse(prob_DMP >= prob & m_diff_tn <= (-effect_size), "hypo", NA))]

  # checkpoint
  print(paste("DMPs annotated given prob = ", prob, " and minimum effect-size = ", 
              effect_size, sep=" ")) 
                          
  # save CAMDAC results per CpG
  output_file1 = file.path(path_output, "CAMDAC_results_per_CpG.RData")
  CAMDAC_results_per_CpG <- dt
  save(CAMDAC_results_per_CpG, file=output_file1)
  
  # checkpoint
  cat(paste0("CAMDAC CpG-wise results saved at:\n", output_file1,
             "\nThis includes copy number information, pure tumour ",
             "methylation rates and DMP calls.\n")) 
  rm(output_file1)

  # extract DMPs
  cols1 <- c("chrom", "start", "end", "DMP_t", "prob_DMP", "m_diff_tn")
  cols2 <- c("chrom", "start", "end", "name", "score", "effect_size")
  CAMDAC_DMPs <- CAMDAC_results_per_CpG[!is.na(DMP_t), .SD, .SDcols=cols1]
  colnames(CAMDAC_DMPs) <- cols2
  rm(CAMDAC_results_per_CpG, cols1, cols2)

  # save DMPs in standard BED5 format
  output_file2 = file.path(path_output, "CAMDAC_DMPs.bed")
  write.table(CAMDAC_DMPs, file=output_file2, sep='\t', col.names = TRUE, quote=FALSE)
  
  # checkpoint
  cat(paste0("CAMDAC DMPs saved in BED5 format at ", output_file2, "\n")) 
  rm(CAMDAC_DMPs, output_file2)  

  # Extract DMPs and obtain summary stats
  tmp <- dt[!is.na(prob_DMP) & prob_DMP>prob,]
  DMP_counts_t <- tmp[, sum(abs(m_t_corr-m_n) > 0.2, na.rm=TRUE)]
  if(bulk == TRUE){ DMP_counts_b <- tmp[, sum(abs(m_b-m_n) > 0.2, na.rm=TRUE)]}
  rm(tmp)

  # set DMP stats var names
  nam <- c("DMP_counts_", "DMP_prop_")
  nams <- paste0(nam, "t")
  vec <- c(DMP_counts_t, round2(DMP_counts_t/nrow(dt), digits=3))
  
  if(bulk == TRUE){
    nams <- c(paste0(nam, "b"), nams)
    vec <- c(DMP_counts_b, round2(DMP_counts_b/nrow(dt), digits=3), vec)
  } 

  # save DMP stats file
  tmp <- data.frame(y=vec, stringsAsFactors=FALSE)
  rownames(tmp) <- nams
  output_file3 <- paste0(path_output,"/", patient_id, "_", sample_id, "_DMP_stats.txt")
  write.table(tmp , file=output_file3, sep='\t', col.names = FALSE, quote=FALSE)

  # checkpoint
  cat(paste0("\nDMP summary stats saved in ",output_file3,"\n"))
  rm(nam, nams, vec, tmp, output_file3)
 
  # load annotations
  annotations_file = paste0(path_to_CAMDAC, "/pipeline_files/",
                            build, "_annotations/", build, "_all_regions_annotations.fst")
  all_regions_anno <- fst::read_fst(path = annotations_file, as.data.table = TRUE)                         
  all_regions_anno[, chrom := factor(chrom, levels=c(1:22,"X","Y"), ordered=TRUE)]
  
  ## Group CpGs into bins, get bin methylation info and annotate Ensembl features
  #dt_bins <- bin_CpGs(path=path_output, patient_id=patient_id, sample_id=sample_id,
  #                   dt=dt, anno_list=all_regions_anno, n_cores=n_cores)
  #colnames(dt_bin)[grepl("^i\\.", colnames(dt_bin))] <- c("start", "end")
  #save(dt_bins, file =paste0(path,"/", "CAMDAC_results_per_bin.RData"))

  # Should both bulk and pure tumour files be generated?
  if(bulk == FALSE){
    prefixes = "purified_tumour"
    flag = FALSE
  } else {
    prefixes = c("purified_tumour", "bulk_tumour")
    flag = c(FALSE, TRUE)
  }

  for(bulk in flag){
      # get DMRs
      CAMDAC_DMRs <- get_DMRs(path=path_output, patient_id=patient_id, sample_id=sample_id, dt=dt,
                              anno_list=all_regions_anno, min_DMP_counts=min_DMP_counts_in_DMR, 
                              min_consec_DMP=min_consec_DMP_in_DMR,
                              n_cores=n_cores, bulk=bulk)
      if (is.null(CAMDAC_DMRs)){
        return(NULL)
      }
      colnames(CAMDAC_DMRs)[grepl("^i\\.", colnames(CAMDAC_DMRs))] <- c("start", "end")
  
      # set filenames and filepaths
      prefix <- "CAMDAC"
      if(length(prefixes)==2){prefix <- paste(prefix, prefixes[which(flag==bulk)], sep="_") }
      suffix <- "annotated_DMRs.RData"
      f_nm <- paste(prefix, suffix, sep="_")
      output_file4 = paste0(path_output, f_nm)
 
      # save DMRs with annotation
      save(CAMDAC_DMRs, file=output_file4)
      rm(dt, f_nm)

      # checkpoint
      cat(paste0(prefix, " DMRs identified and saved at ", output_file4, "\n")) 
      rm(output_file4)

      # set filenames and filepaths
      suffix <- "DMRs.bed"
      f_nm = paste(prefix, suffix, sep="_")
      output_file5 = paste0(path_output, f_nm)
      
      # make temporary column m_diff
      if(bulk == FALSE){ CAMDAC_DMRs[, "m_diff" := m_t - m_n]}
      if(bulk == TRUE){ CAMDAC_DMRs[, "m_diff" := m_b - m_n]}

      # compare bulk, tumour and normal methylomes at DMPs
      plot_methylation_info_with_anno(dt=CAMDAC_DMRs, path=path_output, bulk=bulk)

      # extract essential columns
      cols1 <- c("chrom", "start", "end", "DMR_type", "prob", "m_diff")
      cols2 <- c("chrom", "start", "end", "name", "score", "effect_size")
      CAMDAC_DMRs <- CAMDAC_DMRs[, .SD, .SDcols=cols1]
      colnames(CAMDAC_DMRs) <- cols2

      # save DMPs in standard BED5 format
      write.table(CAMDAC_DMRs, file=output_file5, sep='\t', col.names = TRUE, quote=FALSE)
      rm(cols1, cols2, output_file5)

  }
  # remove pure tumour methylation rates files
  file.remove(paste0(path_output, "purified_tumour", ".RData"))
  rm(all_regions_anno)
}

#' Get DMPs
#' 
#' \code{get_DMPs} returns a df with annotated statistics for each CpG
#'
#' @param path Complete path to the CAMDAC methylation output directory
#' fir this sample
#' @param patient_id Character string containting the patient number
#' @param sample_id Character variable with the tumour sample_id
#' @param df A data.table with pure, bulk and normal methylation info
#' @param prob Numerical  value representing the threshold for statistically
#' significant DMP (default is p=0.99)
#' @param n_cores Number of cores to do the statistical testing over
#' 
#' @return A data.table object with all the CpG loci, their coverage, counts 
#' methylated and methylation rate
get_DMPs <- function (path, patient_id, sample_id, df, prob=0.99, n_cores) {
  
  # Evan Miller's closed form solution for the probability that
  # a draw from a beta dist is dfeater than anoter,
  # p_A ~ Beta(alpha_A, beta_A), p_B ~ Beta(alpha_B, beta_B).
  # In this case, A is the normal and B is the bulk tumour.
  # Alpha is the counts methylated, beta is counts unmethylated.
  h <- function(alpha_n, beta_n, alpha_b, beta_b){
    j <- seq.int(0, round(alpha_b)-1)
    log_vals <- (lbeta(alpha_n + j, beta_n + beta_b) - log(beta_b + j) -
                   lbeta(1 + j, beta_b) - lbeta(alpha_n, beta_n))
    1 - sum(exp(log_vals))
  }
    
  #Define Variables and add pseudocount
  alpha_n=df$M_n+1
  beta_n=df$UM_n+1
  alpha_b=df$M_b+1
  beta_b=df$UM_b+1
  # Although a B(0.5, 0.5) prior is more suited, to methylation rates, 
  # alpha_b must be an integer value for h() to hold.
  # We therefore add a pseudocount of 1, meaning that B goes to B(1, 1) as
  # alpha and beta values go to zero).

  # set up results object
  n <-length(df$M_b)
  result <- cbind(numeric(n))

  # Get DMPs
  result[,1] <- parallel::mcmapply(function(alpha_n,beta_n,alpha_b,beta_b) {
    prob_hypo <- NULL
    prob_hypo <- h(alpha_n = alpha_n, beta_n=beta_n, alpha_b=alpha_b, beta_b=beta_b)
    if(is.null(prob_hypo)){prob_hypo <- NA}
    return(prob_hypo)
  }, alpha_n=alpha_n,beta_n=beta_n,alpha_b=alpha_b,beta_b=beta_b,mc.cores=n_cores)
  
  # Cap and flip probabilities where necessary
  tmp_m_b <- alpha_b/(alpha_b+beta_b)
  tmp_m_n <- alpha_n/(alpha_n+beta_n)
  df[, prob_DMP := ifelse(result[,1]>1, 1, ifelse(result[,1]<0, abs(result[,1]), result[,1]))]
  df[, prob_DMP := ifelse(tmp_m_b - tmp_m_n > 0, 1-prob_DMP, prob_DMP)]
 
  # return results object
  return(df)
  
  # clean-up
  rm(result,n,alpha_n,beta_n,alpha_b,beta_b)
}


#' Cluster CpGs into annotated bins
#' 
#' \code{bin_CpGs} returns the df with the annotation for each CpG
#' 
#' @param path Character string of the output directory
#' @param patient_id Character string containting the patient ID
#' @param sample_id Character string containting the sample ID. 
#' @param dt data.table where each CG is a row with DMP info. 
#' @param anno_list A data.table object containing annotated genomic bins including
#' genes, exons, introns, UTRs, CGI, CGI shores, CGI shelves, promoters or enhancers
#' @param n_cores number of cores for parallel processing
#'  
#' @return A dataframe for each sample_id with the copy number calls added
bin_CpGs <- function (path, patient_id, sample_id, dt, anno_list, n_cores) {

  # Ensure dt is a data.table object
  if(class(dt)[1]!="data.table"){
  dt <- data.table(dt)
  }
  
  # ensure both seqlevels to characters
  dt[, chrom := as.character(chrom)]
  anno_list[, chrom := as.character(chrom)]
  
  # overlap annotated regions and CpG methylation objects
  setkey(dt, chrom, start, end)
  setkey(anno_list, chrom, start, end)
  ov <- foverlaps(anno_list, dt, type="any", nomatch=NULL)

  # order the test regions
  ov <- ov[order(as.numeric(cluster_id),
                 factor(chrom, levels=c(1:22, "X", "Y"), ordered=TRUE),
                 start,end),]
  
  # extract all bin ids with coverage and number of unique bins
  ids <- unique(ov$cluster_id)
  l = length(ids)
 
  cat("Concatenate annotated bins\n")
  # concatenate annotated CpG methylation
  dt_anno_bins <- rbindlist(parallel::mclapply(1:l, function(i, df, ids){
      x <- df[cluster_id==ids[i], ]
      x <- x[, segment := paste0(chrom,":",seg_start,"-",seg_end)]
      y <- x[, .(m_n= mean(m_n, na.rm=TRUE),
                 cov_n= mean(cov_n, na.rm=TRUE), 
                 m_n_low= mean(m_n_low, na.rm=TRUE),  
                 m_n_high= mean(m_n_high, na.rm=TRUE),
                 m_b= mean(m_b, na.rm=TRUE),
                 cov_b= mean(cov_b, na.rm=TRUE),  
                 m_b_low= mean(m_b_low, na.rm=TRUE),  
                 m_b_high= mean(m_b_high, na.rm=TRUE), 
                 m_t= mean(m_t_corr, na.rm=TRUE),
                 cov_t= mean(cov_t, na.rm=TRUE), 
                 m_t_low= mean(m_t_low, na.rm=TRUE),  
                 m_t_high= mean(m_t_high, na.rm=TRUE), 
                 prob= mean(prob_DMP, na.rm=TRUE),
                 CN= mean(CN, na.rm=TRUE), 
                 nA= mean(nA, na.rm=TRUE), 
                 nB= mean(nB, na.rm=TRUE),
                 segment= paste(unique(segment), collapse = ";")),
             by = .(cluster_id, chrom, i.start, i.end)]
      return(y)
      }, df=ov, ids=ids, mc.cores=n_cores))

  # Add bin annotations
  colls <- c("intragenic", "gene_ids", "CGI", "enhancer", "enhancer_ids",
             "promoter", "gene_prom_ids", "exon", "intron", "FUTR", "TUTR", 
             "repeats", "repeat_name", "repeat_class", "repeat_family")
  dt_anno_bins[, as.character(colls) :=
                 anno_list[match(dt_anno_bins$cluster_id, cluster_id), .SD, .SDcols=colls]]

  # return
  return(dt_anno_bins)
} 

#' Assign bins 
#' 
#' \code{annotate_DMRs} returns the df with the annotation for each CpG
#' 
#' @param path Character string of the output directory
#' @param patient_id Character string containting the patient_id ID
#' @param sample_id Character string containting the sample ID. 
#' @param dt dataframe where each CG is a row with DMP info. 
#' @param anno_list A data.table object containing annotated genomic bins including
#' genes, exons, introns, UTRs, CGI, CGI shores, CGI shelves, promoters or enhancers
#' @param min_DMP_counts Numerical - number of DMPs required in a DMR
#' @param min_consec_DMP Numerical - number of consecutive DMPs required in a DMR
#' @param n_cores number of cores for parallel processing
#'  
#' @return A dataframe for each sample_id with the copy number calls added
get_DMRs <- function (path, patient_id, sample_id, dt, anno_list,
                      min_DMP_counts, min_consec_DMP, n_cores, bulk=FALSE) {

  # ensure both seqlevels to characters
  dt[, chrom := as.character(chrom)]
  anno_list[, chrom := as.character(chrom)]
  
  # overlap annotated regions and CpG methylation objects
  setkey(dt, chrom, start, end)
  setkey(anno_list, chrom, start, end)
  ov <- foverlaps(anno_list, dt, type="any", nomatch=NULL)

  # order the test regions
  ov <- ov[order(as.numeric(cluster_id),
                 factor(chrom, levels=c(1:22, "X", "Y"), ordered=TRUE),
                 start,end),]
   
  # print analysis parameters
  cat(paste0("DMR threholds","\n", "min DMP counts : ", min_DMP_counts ,"\n", 
             "min number of consecutive DMPs : ", min_consec_DMP, "\n"))  
  
  # annotations to be assigned
  anno_names <- c("all_CpGs", "intergenic", "intragenic", "CGI", "shore", "shelf", 
                  "enhancer", "promoter", "exon", "intron", "FUTR", "TUTR", "repeats")

  # store DMR stats output file              
  output_file <- paste0(path, patient_id, "_", sample_id, "_DMR_stats.txt")

  # Create data.frame to store DMR stats 
  if(bulk==FALSE){
  DMR_stats <- data.frame(num_CpGs_annot = numeric(length=length(anno_names)), 
                          num_bins = 0, mean_CpGs_per_bin = 0, mean_m_n_DMRs = 0,
                          frac_DMRs_t = 0, mean_m_t_DMRs = 0,                                  
                        stringsAsFactors = FALSE)
  rownames(DMR_stats) <- anno_names
  } else {
    DMR_stats <- as.data.frame(read.delim())
    DMR_stats <- cbind(DMR_stats, frac_DMRs_b = 0,mean_m_b_DMRs = 0) 
  }  

  # get total number of individual CpGs and CpG cluters
  if(bulk==FALSE){
  DMR_stats$num_CpGs_annot[1] <- nrow(dt)
  DMR_stats$num_bins[1] <- length(unique(ov$cluster_id))
  }

  # find consectuvive DMPs
  fct <- function(x,sgn,min) {
    flag <- c(FALSE, diff(sign(sgn)) > 0)
    y <- rle(!is.na(x) & flag==FALSE)
    z<-as.numeric(suppressWarnings(
      max(y$lengths[y$lengths >= min & y$values == TRUE]))) 
    if(!is.finite(z)){z <- as.numeric(NA)}
    return(z)
  } 

  # note relevant columns names
  if(bulk==FALSE){
    colls <- c("m_t_corr","m_n","DMP_t","cluster_id")
  } else {
    colls <- c("m_b","m_n","DMP_b","cluster_id")
  }
  colls2 <- c("m","m_n","DMP","cluster_id")

  # extract and re-name columns
  tmp <- ov[, .SD, .SDcols=colls]
  colnames(tmp) <- colls2

  # Get CpG and DMP counts for each cluster
  results <- tmp[,.(CpG_counts = length(m), 
                    DMP_counts = length(m[!is.na(DMP)])),
                   by=.(cluster_id)][order(as.numeric(cluster_id)),] 
    
  # Get counts of consecutive DMPs for each cluster
  tmp2 <- tmp[,.(consec_DMP = fct(x = DMP, sgn=m-m_n, min = min_consec_DMP)),
                by=.(cluster_id)] 
  results$consec_DMP <- tmp2[match(as.numeric(results$cluster_id),as.numeric(cluster_id)), 
                             .SD, .SDcols=c("consec_DMP")] 
  rm(tmp, tmp2, fct)
  
  # Annotate DMR calls depending on consecutive DMP and min DMP counts tresholds
  results[, DMR := ifelse(DMP_counts >= min_DMP_counts & 
                          consec_DMP >= min_consec_DMP,"DMR",NA)]
  
  # store total DMRs stats
  if(bulk==TRUE){ k <- which(colnames(DMR_stats)=="frac_DMRs_b") }
  if(bulk==FALSE){ 
    k <- which(colnames(DMR_stats)=="frac_DMRs_t") 
    DMR_stats$mean_CpGs_per_bin[1] <- round2(results[,mean(CpG_counts)], digits = 0)
  }
  DMR_stats[1, k] <- round2(sum(!is.na(results$DMR))/nrow(results), digits=3)
  
  # extract all bin ids with coverage and number of unique bins
  ids <- results[!is.na(DMR), unique(cluster_id)]
  l <- length(ids)

  # Report and return if no DMRs found
  if (length(ids) == 0){
    cat("No DMRs found with the current parameters.\n")
    return(NULL)
  } else {
    cat(paste0("Number of DMRs found: ", length(ids), "\n"))
  }
 
  cat("Concatenate DMR calls \n")
  # concatenate annotated CpG methylation
  if(bulk==FALSE){
    dt_DMRs <- rbindlist(parallel::mclapply(1:l, function(i, df, ids){
      x <- df[cluster_id==ids[i], ]
      x <- x[, segment := paste0(chrom,":",seg_start,"-",seg_end)]
      y <- x[, .(m_n= ifelse(sum(!is.na(DMP_t))==0, as.numeric(NA), mean(m_n[!is.na(DMP_t)])),
                 m_n_low= ifelse(sum(!is.na(DMP_t))==0, as.numeric(NA), mean(m_n_low[!is.na(DMP_t)])),
                 m_n_high= ifelse(sum(!is.na(DMP_t))==0, as.numeric(NA), mean(m_n_high[!is.na(DMP_t)])),
                 m_t= ifelse(sum(!is.na(DMP_t))==0, as.numeric(NA), mean(m_t_corr[!is.na(DMP_t)])),
                 m_t_low= ifelse(sum(!is.na(DMP_t))==0, as.numeric(NA), mean(m_t_low[!is.na(DMP_t)])),
                 m_t_high= ifelse(sum(!is.na(DMP_t))==0, as.numeric(NA), mean(m_t_high[!is.na(DMP_t)])),
                 prob= ifelse(sum(!is.na(DMP_t))==0, as.numeric(NA), mean(prob_DMP[!is.na(DMP_t)])),
                 CG_CN= mean(CG_CN, na.rm=TRUE), nA= mean(nA, na.rm=TRUE), nB= mean(nB, na.rm=TRUE),
                 segment= paste(unique(segment), collapse = ";"),
                 DMR_type= ifelse(sum(!is.na(DMP_t))==0, "NA",
                           ifelse(sum(is.na(DMP_t) | DMP_t=="hypo") >= 0.9*length(m_n), "hypo", 
                           ifelse(sum(is.na(DMP_t) | DMP_t=="hyper") >= 0.9*length(m_n), "hyper", 
                                  "mixed")))),
             by = .(cluster_id, chrom, i.start, i.end)]
      return(y)
    }, df=ov, ids=ids, mc.cores=n_cores))
  } else {
    dt_DMRs <- rbindlist(parallel::mclapply(1:l, function(i, df, ids){
      x <- df[cluster_id==ids[i], ]
      x <- x[, segment := paste0(chrom,":",seg_start,"-",seg_end)]
      y <- x[, .(m_b= ifelse(sum(!is.na(DMP_b))==0, as.numeric(NA), mean(m_b[!is.na(DMP_b)])),
                 m_b_low= ifelse(sum(!is.na(DMP_b))==0, as.numeric(NA), mean(m_b_low[!is.na(DMP_b)])),
                 m_b_high= ifelse(sum(!is.na(DMP_b))==0, as.numeric(NA), mean(m_b_high[!is.na(DMP_b)])),
                 m_n=  ifelse(sum(!is.na(DMP_b))==0, as.numeric(NA), mean(m_n[!is.na(DMP_b)])),
                 m_n_low= ifelse(sum(!is.na(DMP_b))==0, as.numeric(NA), mean(m_n_low[!is.na(DMP_b)])),
                 m_n_high= ifelse(sum(!is.na(DMP_b))==0, as.numeric(NA), mean(m_n_high[!is.na(DMP_b)])),
                 prob= ifelse(sum(!is.na(DMP_b))==0, as.numeric(NA), mean(prob_DMP[!is.na(DMP_b)])),
                 CG_CN= mean(CG_CN, na.rm=TRUE), nA= mean(nA, na.rm=TRUE), nB= mean(nB, na.rm=TRUE),
                 segment= paste(unique(segment), collapse = ";"),
                 DMR_type= ifelse(sum(!is.na(DMP_b))==0, "NA",
                           ifelse(sum(is.na(DMP_b) | DMP_b=="hypo") >= 0.9*length(m_n), "hypo", 
                           ifelse(sum(is.na(DMP_b) | DMP_b=="hyper") >= 0.9*length(m_n), "hyper",
                                  "mixed")))),
             by = .(cluster_id, chrom, i.start, i.end)]
      return(y)
    }, df=ov, ids=ids, mc.cores=n_cores))
  }
  rm(l, ids)

  # Add DMRs annotation and related info to the annotated bins
  dt_DMRs[, c("CpG_counts", "DMP_counts", "consec_DMPs", "DMR") := 
            c(results[match(dt_DMRs$cluster_id, cluster_id), 
                .(CpG_counts, DMP_counts, consec_DMP, DMR)])]

  # store methylation rates stats
  if(bulk==TRUE){   
  DMR_stats$mean_m_b_DMRs[1] <-  # subset to DMP hotspots in DMRs
        round2(dt_DMRs[!is.na(DMR), mean(m_b)], digits = 3)
   } else {
  DMR_stats$mean_m_n_DMRs[1] <-  # subset to DMP hotspots in DMRs
        round2(dt_DMRs[!is.na(DMR), mean(m_n)], digits = 3)
  DMR_stats$mean_m_t_DMRs[1] <-  # subset to DMP hotspots in DMRs
        round2(dt_DMRs[!is.na(DMR), mean(m_t)], digits = 3)
    }

  # Add bin annotations
  colls <- c("intragenic", "gene_ids", "CGI", "enhancer", "enhancer_ids",
             "promoter", "gene_prom_ids", "exon", "intron", "FUTR", "TUTR", 
             "repeats", "repeat_name", "repeat_class", "repeat_family")
  dt_DMRs[, as.character(colls) :=
                 anno_list[match(dt_DMRs$cluster_id, cluster_id), .SD, .SDcols=colls]]
  
  for(anno in anno_names[2:length(anno_names)]){
    # extrat row id
    j = which(anno_names==anno)
    k <- ifelse(anno=="intergenic", which(colnames(dt_DMRs)=="intragenic"), 
         ifelse(anno%in%c("shore", "shelf"), which(colnames(dt_DMRs)=="CGI"), 
           which(colnames(dt_DMRs)==anno)))

    # extrac the relevant list
    vec <- unname(unlist(dt_DMRs[, .SD, .SDcols=c(k)]))
    flag <- vec==anno
    rm(vec)
    
    # get total number of CpGs, DMPs and regions with coverage
    DMR_stats$num_CpGs_annot[j] <- dt_DMRs[flag==TRUE, sum(CpG_counts, na.rm=TRUE)]
    DMR_stats$num_bins[j] <- dt_DMRs[flag==TRUE, length(unique(cluster_id))]
   
    # store bin stats
    DMR_stats$mean_CpGs_per_bin[j] <- 
        round2(dt_DMRs[flag, mean(CpG_counts)], digits = 0)
 
    # store DMR stats
    if(bulk==TRUE){ k <- which(colnames(DMR_stats)=="frac_DMRs_b") }
    if(bulk==FALSE){ k <- which(colnames(DMR_stats)=="frac_DMRs_t") }
    DMR_stats[j, k] <- dt_DMRs[flag, sum(!is.na(DMR))] / nrow(dt_DMRs[flag,])
    
    # store bin/DMR mean methylation levels
    DMR_stats$mean_m_n_DMRs[j] <-  # subset to DMP hotspots
        round2(dt_DMRs[flag & !is.na(DMR), mean(m_n)], digits = 3)
    if(bulk==TRUE){   
        DMR_stats$mean_m_b_DMRs[j] <-  # subset to DMP hotspots
            round2(dt_DMRs[flag & !is.na(DMR), mean(m_b)], digits = 3)
    } else {
        DMR_stats$mean_m_t_DMRs[j] <-  # subset to DMP hotspots
            round2(dt_DMRs[flag & !is.na(DMR), mean(m_t)], digits = 3)
    }
  }
  rm(anno, anno_names, j, k)
  
  # save stats
  write.table(DMR_stats, file=output_file, sep='\t', col.names = TRUE, quote=FALSE) 
  rm(DMR_stats, output_file)

  # return
  return(dt_DMRs)
} 

# Plot summary methylation information with annotated information
# Arguments:
#' @title Plot methylation information
#' @param dt Data table with methylation information per CpG
#' @param path Character path variable pointing to the desired working directory.
#' @param bulk Logical determining whether the bulk or purified tumour is to be plotted
#' @return NULL
#' @keywords internal
plot_methylation_info_with_anno <- function(dt, path, bulk){

 # Set color code for hyper/hypo 
 cols <- c("gold2", "royalblue2")
 names(cols) <- c("hyper", "hypo")
 
 # Set color code for annotations
 cols2 <- c("black", "grey40", "grey75", "red2", "pink2", "dodgerblue", "forestgreen", 
            "olivedrab3", "lightgreen", "mediumpurple")
 names(cols2) <- c("all_bins", "exon", "intron", "promoter", "prom_CGI", "enhancer","CGI", 
                   "shore", "shelf","repeats")

  tmp <- dt
  if(bulk==FALSE){
       prefix = "purified"
      ylabb = "pure tumour - normal methylation rate"
   } else {
      prefix = "bulk"
      ylabb = "bulk tumour - normal methylation rate"
  }

 # plot fraction hyper/hypo DMRs
 p <- ggplot(data.frame(cols2), aes(x=names(cols2))) + theme_classic() +
      scale_x_discrete(name="", breaks = names(cols2), limits=names(cols2),
                   labels=c("all\nCpG bins", "exon", "intron", "promoter", 
                            "prom\nCGI", "enhancer", "CGI", "shore", "shelf", "repeat")) +
      theme(axis.text.x = element_text(angle = 40, size = 10, vjust=0.75), axis.title.x = element_blank(),
            axis.title.y = element_text(size = 10))
 tmp2 <- data.table(x=names(cols2), y1=0, y2=0, y3=0)
 for(anno in names(cols2)){
   if(anno=="all_bins"){
      flag <- rep(TRUE, nrow(tmp))
   } else if(anno=="prom_CGI"){ 
      flag <- tmp[, promoter=="promoter"&CGI=="CGI"]
   } else {
      nam <- ifelse(anno=="shore", "CGI", 
             ifelse(anno=="shelf", "CGI", 
             ifelse(anno=="repeats", "repeats", anno)))  
       flag <- as.logical(tmp[, .SD, .SDcols=nam]==anno)
   }
   z1 <- tmp[flag, sum(!is.na(DMR_type)&DMR_type=="hyper")/sum(!is.na(DMR_type))]
   z2 <- tmp[flag, sum(!is.na(DMR_type)&DMR_type=="hypo")/sum(!is.na(DMR_type))]
   tmp2[x==anno, c("y1", "y2") := list(z1, z2)]
 
 }
 p1 <- p + scale_y_continuous(name= "hypo-/hypermethylated DMR ratio", breaks = seq(-1,1,0.2), limits=c(-1, 1)) +
           geom_col(data=tmp2, aes(x=x,y=y1,col="hyper", fill="hyper"), alpha=0.5) +
           geom_col(data=tmp2, aes(x=x,y=(-1+y1),col="hypo", fill="hypo"), alpha=0.5) +
           scale_fill_manual(name="anno",values=cols) + 
           scale_color_manual(name="anno",values=cols)

 # plot methylation rate difference per anno
 p2 <- p + geom_hline(yintercept=0, size=0.25) + scale_y_continuous(name= ylabb, breaks = seq(-1,1,0.2)) +
  geom_violin(data=tmp, aes(x="all_bins", y=m_diff, col="all_bins", fill="all_bins"), alpha=0.2) +
  geom_violin(data=tmp[exon=="exon",], aes(x="exon",y=m_diff, col="exon", fill="exon"), alpha=0.2) +
  geom_violin(data=tmp[intron=="intron",], aes(x="intron",y=m_diff, col="intron", fill="intron"), alpha=0.2) +
  geom_violin(data=tmp[promoter=="promoter",], aes(x="promoter",y=m_diff, col="promoter", fill="promoter"), alpha=0.2) +
  geom_violin(data=tmp[promoter=="promoter"&CGI=="CGI",], aes(x="prom_CGI",y=m_diff, col="prom_CGI", fill="prom_CGI"), alpha=0.2) +
  geom_violin(data=tmp[enhancer=="enhancer",], aes(x="enhancer",y=m_diff, col="enhancer", fill="enhancer"), alpha=0.2) +
  geom_violin(data=tmp[CGI=="CGI",], aes(x="CGI", y=m_diff, col="CGI", fill="CGI"), alpha=0.2) +
  geom_violin(data=tmp[CGI=="shelf",], aes(x="shelf", y=m_diff, col="shelf", fill="shelf"), alpha=0.2) +
  geom_violin(data=tmp[CGI=="shore",], aes(x="shore", y=m_diff, col="shore", fill="shore"), alpha=0.2) +
  geom_violin(data=tmp[repeats=="repeats",], aes(x="repeats",y=m_diff, col="repeats", fill="repeats"), alpha=0.2) +
  #geom_jitter(width=0.25, size=2, shape=20, alpha=0.9) +  
  scale_fill_manual(name="anno",values=cols2) + 
  scale_color_manual(name="anno",values=cols2) 
  
  gglist1 <- list(p1)
  gglist2 <- list(p2)
  my_layout <- rbind(1, 2)
  gd <- arrangeGrob(grobs = c(gglist1, gglist2), layout_matrix = my_layout)
  Sys.sleep(1)
  ggsave(gd, filename = paste0(path, "/", prefix, "_tumour_DMR_summary_plots.pdf"), device="pdf",
             width=9, height = 6, unit="in")
}

#' @description Round numerical values to 'n' digits
#' @param x Numerical vector containing the numbers to round
#' @param digits Numerical value representing the number of decimal digits to retain
#' @return rounded numerical vector
round2 <- function(x, digits) {
  ifelse(as.integer(x*(10^(digits+1))) %% 10 >= 5,
         ceiling(x*(10^digits))/(10^digits),
         floor(x*(10^digits))/(10^digits))
}

# END
