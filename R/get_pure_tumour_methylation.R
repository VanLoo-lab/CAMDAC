##  Get purified tumour methylation rates from bulk tumour
##  bisulfite sequencing data
##  Latest update: 04/05/2021
##  Version 1.0.0
##  Author: Elizabeth Larose Cadieux

#' Deconvolve the pure tumour methylation rate from bisulfite sequencing data
#'
#' \code{get_pure_tumour_methylation}
#' 
#' @param patient_id Character variable containting the patient id number
#' @param sample_id Character variable with the (control or tumour) sample_id
#' @param sex Character variable with  the patient expressed as "XX" for
#' female or "XY" for male.
#' @param normal_infiltrates_proxy_id, Sample ID of the matched normal control
#' @param path Character path variable pointing to the desired working directory.
#' This is where the output will be stored and should be constant for all CAMDAC functions.
#' @param path_to_CAMDAC Character variable containting the path to the CAMDAC
#' directory including dir name (e.g. "/path/to/CAMDAC/").
#' @param build Character variable corresponding to the reference genome
#' used for alignment. CAMDAC is compatible with "hg19", "hg38", "GRCH37","GRCH38".
#' @param n_cores Numerical value correspdonding to the number of cores
#' for parallel processing
#' @param reseg Logical value should be set to FALSE. Multi-sample re-segmentation of
#' the copy number profiles will be available in future versions of CAMDAC.
#'
#' Note:
#' #' Annotation include:
#' CGI (including shore and shelves)
#' gene body (intragenic, 5UTR, 3UTR, intron, exon)
#' promoter (2kb upstream and 500 downstream any UCSC annotated gene)
#' enhancer (vista and FANTOM5 annotation)
#'
#' @return CAMDAC purified tumour methylation rates

get_pure_tumour_methylation <- function(patient_id,sample_id,sex,
                                        normal_infiltrates_proxy_id,
                                        path,path_to_CAMDAC,build,
                                        n_cores,reseg=FALSE){

  control_id=normal_infiltrates_proxy_id
  if(sample_id==control_id){
      stop(paste("Purified tumour methylation rate cannot be obtained for control",
                 "samples as this sample is a proxy for the normal methylation rate",
                 sep="\n"))
  }
  if(detectCores()<n_cores){
      warning(paste0(n_cores, " cores selected but only ",
      detectCores(), " detected on machine."))
      
  }

  # set build to to assembly version disregarging UCSC/Ensembl
  if(build=="GRCH37"){build="hg19"}
  if(build=="GRCH38"){build="hg38"}
  
  # set output path
  path_patient <- paste0(file.path(path, patient_id),"/")
  path_output <-paste0(path_patient, "Methylation/",sample_id,"/")
  if(reseg == TRUE){ path_output <- paste0(path_output, "/re_segmented/") }
    
  # set filenames and filepaths
  # annotations_file = paste0(path_to_CAMDAC, "pipeline_files/",
  # build, "_annotations/", build, "_all_regions_annotations.RData")
    
  # Extract methylation data
  file_methylation = paste0(path_patient, "/Methylation/", sample_id, "/dt_tumour_and_normal_m.RData")
  load(file_methylation)
  dt_mb_mn <- dt_tumour_and_normal_m
  rm(dt_tumour_and_normal_m)

  # Load germline genotypes
  germline_genotype_filename = paste0(path_patient,"/Copy_number/",sample_id,"/",
                                      patient_id, ".",sample_id,".SNPs.RData")
  load(germline_genotype_filename)
  rm(germline_genotype_filename)
  
  # Get heterosygous SNPs
  dt_SNPs <- dt[type=="Heterozygous",] 
  dt_SNPs[, c("CHR", "start", "end") := .(paste0("chr", chrom), POS, POS)]
  setkeyv(dt_SNPs, cols=c("CHR", "start", "end"))
  rm(dt)
  
  # Find heterozygous SNPs overlapping with CpGs 
  dt_mb_mn[, offset := ifelse((end-start)>2,1,0)] 
  dt_mb_mn[, ID := 1:nrow(dt_mb_mn)]
  tmp <- dt_mb_mn[, .(ID, CHR, start=start+offset, end=end-offset)]
  setkeyv(tmp, cols=c("CHR", "start", "end"))
  ov <- foverlaps(dt_SNPs, tmp, nomatch = 0)

  # Label CG-forming and CG-destroying SNPs
  ov[, c("i.start", "i.end") := NULL]
  ov[, label := ifelse((POS==start & ref=="C")|(POS==end & ref=="G"),
                       "CG-destroying","CG-forming")]
  
  # store label and tumour BAF for downstream purified tumour methylation rate calculation
  dt_mb_mn[, c("BAF","label") := ov[match(dt_mb_mn$ID,ov$ID), .(BAF,label)]]
  dt_mb_mn[, c("offset", "ID") := NULL]  
  rm(dt_SNPs, tmp)
  
  # Load the copy number segments file
  if(reseg == FALSE){
    file_CN = paste0(path_patient, "/Copy_number/", sample_id, "/raw/", patient_id,
                     ".",sample_id,".ascat.output.RData")
        if(!file.exists(file_CN)){
          file_CN = paste0(path_patient, "/Copy_number/", sample_id, "/", patient_id,
                         ".",sample_id,".ascat.output.RData")
        }
    }
  if(reseg == TRUE){
    file_CN = paste0(path_patient, "/Copy_number/", sample_id, "/re_segmented/",
                     patient_id, ".",sample_id,".ascat.output.resegmented.RData")
    }
  load(file_CN)
    
  # Get copy number segments and tumour purity
  seg <- ascat.output$segments_raw
  p = ascat.output$aberrantcellfraction
  rm(ascat.output)
  
  # set chromosome labels to Ensembl style
  colnames(dt_mb_mn)[1] <- "chrom"
  dt_mb_mn[, chrom := substr(chrom, 4, 5)]

  # Run copy number data extraction
  dt_mb_mn_cn <- annotate_copy_number(dt_sample=dt_mb_mn, 
                                      seg=seg,
                                      rm_sex_chrom=FALSE)
  rm(seg, dt_mb_mn)

  # Checkpoint
  print(paste("Copy number estimates added for each CpG locus", sep=" "))
  
  # Annotate allele-specific copy number
  #dt_mb_mn_cn[, multi := paste(nA, nB, sep="+")]
  
  # Tumour purity estimate for each sample
  print(paste0("Sample purity is estimated at rho = ", p,
               ". Copy number segments assigned to CpG loci."))
  
  # Compute purified tumour methylation rates
  dt_mt <- compute_tumour_methylome(dt = dt_mb_mn_cn, 
                                    p = p, min_cov_t = 3,
                                    sex=sex, build=build)
  print("Purified tumour methylation rate estimates added for each CpG locus")
  rm(dt_mb_mn_cn)

  # Set vars for m_t HDI calculation 
  n <- nrow(dt_mt)
  HDI <- cbind.data.frame(numeric(n), numeric(n))
  M_b=dt_mt$M_b ; UM_b=dt_mt$UM_b
  if("M_n" %in% colnames(dt_mt)){
      M_n=dt_mt$M_n ; UM_n=dt_mt$UM_n
  }
  if("M_n_i" %in% colnames(dt_mt)){
      M_n=dt_mt$M_n_i ; UM_n=dt_mt$UM_n_i
  }
  CN=dt_mt$CG_CN
  CN_n=dt_mt$CG_CN_n

  # Compute the tumour methylation rate HDI 
  HDI[,1:2] <- t(mcmapply(function(M_b,UM_b,M_n,UM_n,p,CN, CN_n)
                          HDIofMCMC(M_b,UM_b,M_n,UM_n,p,CN,CN_n, 0.99),
                          M_b,UM_b,M_n,UM_n,p,CN,CN_n,mc.cores=n_cores))             
  dt_mt[,  c("m_t_low", "m_t_high") := as.list(HDI)]
  dt_mt[, "CG_CN_n" := NULL]
  
  # set output filename
  output_file = paste0(path_output, "purified_tumour", ".RData")
  rm(M_b, UM_b, M_n, UM_n, HDI, CN, CN_n, n)
  
  # Checkpoint
  print("Purified tumour methylation rate HDI added for each CpG locus")

  # Save CAMDAC outputs so far
  dt_purified_tumour <- dt_mt ; rm(dt_mt)
  save(dt_purified_tumour, file=output_file)

  # Subset to essential columns
  cols1 <- c("chrom", "start", "end","m_t_corr")
  cols2 <- c("chrom", "start", "end", "score")
  tmp <- dt_purified_tumour[, .SD, .SDcols=cols1]
  colnames(tmp) <- cols2

  # save DMPs in standard BED4 format
  output_file = file.path(path_output, "CAMDAC_purified_tumour.bed")
  write.table(tmp, file=output_file, sep='\t', col.names = TRUE, quote=FALSE)

  # checkpoint
  print(paste0("CAMDAC pure tumour methylation rates saved in BED4 format at ", output_file)) 
  rm(tmp, cols1, cols2, output_file)  

  # compare bulk, tumour and normal methylomes
  plot_2d_density(dt=dt_purified_tumour , path=path_output)

  # clean-up redundant files
  file.remove(paste0(path_output, "dt_tumour_and_normal_m.RData"))
}

#' Assign copy number calls
#'
#' \code{annotate_copy_number} returns the data.table dt_sample annotated with allele-specific copy numbers
#'
#' @param dt_sample data.table object with each CpG and their coverage, counts methylated and methylation rate
#' @param seg  ASCAT.m copy number segements object 
#' @param rm_sex_chrom Logical indicating if you would like to remove sex chrom from downstream analyses
#'
#' @return A dataframe for each sample_id with the copy number calls added
annotate_copy_number <- function (dt_sample, seg, rm_sex_chrom=FALSE) {
  
  # Choose the columns that you need to simplify the objects in the subsequent computations
  seg = data.table(chrom = as.character(seg$chr),
                   start = as.numeric(as.character(seg$startpos)),
                   end = as.numeric(as.character(seg$endpos)),
                   nA = seg$nMajor, nB = seg$nMinor,
                   CN = seg$nMajor + seg$nMinor,
                   seg_min = as.numeric(as.character(seg$startpos)),
                   seg_max = as.numeric(as.character(seg$endpos)))
  setkeyv(seg, cols=c("chrom", "start", "end"))
  if(rm_sex_chrom == TRUE){seg <- seg[as.character(seg$chrom)%in%as.character(c(1:22)),]}
  
  # Overlap CpG-wise methylation data.table with copy number segments
  setkeyv(dt_sample, cols=c("chrom", "start", "end"))
  overlaps <- foverlaps(seg, dt_sample, nomatch = 0)

  # Clean up annotated data
  cols <- colnames(dt_sample)
  dt_sample_cn <- cbind(overlaps[, .SD, .SDcols = cols],
                        overlaps[, .(nA, nB, CN, "seg_start"=seg_min, "seg_end"=seg_max)])
  rm(cols)

  # Remove homozygous deletion (CN == 0)
  dt_sample_cn[, FLAG := CN == 0]
  dt_sample_cn <- dt_sample_cn[FLAG == FALSE,]
  dt_sample_cn[, "FLAG" := NULL]

  return(dt_sample_cn)
}


#' Compute the tumour methylation rate
#'
#' \code{compute_tumour_methylome} returns the data.table dt annotated with 
#' CAMDAC pure tumour methylation rates
#'
#' @param dt data.table object with each CpG and their coverage, counts methylated,
#' methylation rate and copy number and matched normal methylation info
#' @param p Numerical - Sample purity estimates
#' @param min_cov_t Numerical - Minimum tumour coverage
#' @param sex Character variable with  the patient expressed as "XX" for female or "XY" for male.
#' @param build Character variable corresponding to the reference genome used for alignment.
#'
#' @return A dataframe for each sample_id with the tumour methylome added
compute_tumour_methylome <- function (dt, p, min_cov_t = 3, sex, build) {

# convert factors to characters
dt[, CN := as.numeric(as.character(CN))]

# Get CG allele total copy number
dt[, CG_CN := ifelse(is.na(BAF), CN,
                  ifelse((BAF<0.5 & label=="CG-destroying")|(BAF>0.5 & label=="CG-forming"), nA,
                         nB))]

# Remove loci with CG allele CN = 0
# This can occur at heterozygous SNPs in regions of LOH or
# at homozygous deletions
dt <- dt[CG_CN>0, ] 

# Get normal copy number at heterozygous SNPs
CN_norm <- ifelse(!is.na(dt$BAF), 1, 2) 
# CN norm is 1 out of 2 at CG-destroying/-forming heterozygous SNPs

# Get normal copy number on sex chromosome X. 
# Chromosome X in males has CN = 1 or 2 respectively outside or within PAR regions.
if(sex=="XY"){
  if(build=="hg19"){
  flag <- dt[, ifelse(chrom=="X" & start>60001 & end<2699520, TRUE, # PAR1
               ifelse(chrom=="Y" & start>10001 & end<2649520, TRUE, # PAR1
               ifelse(chrom=="Y" & start>59034050 & end<59363566, TRUE, # PAR2
               ifelse(chrom=="X" & start>154931044 & end<155260560, TRUE, # PAR2
                      FALSE))))]
  }
  if(build=="hg38"){
  flag <- dt[, ifelse((chrom=="X" | chrom=="Y")& start>10001 & end<2781479,TRUE,
               ifelse(chrom=="Y" & start>56887903 & end<57217415, TRUE,
               ifelse(chrom=="X" & start>155701383 & end<156030895, TRUE,
                      FALSE)))]
    
  }
  CN_norm[flag==FALSE & dt$chrom=="X"] <- 1
  rm(flag)
}

# Get normal CG allele total copy number
dt[, CG_CN_n := CN_norm]

if(sum(grepl("^M_n$", colnames(dt)))==1){
    m_n <- dt$m_n
}
if(sum(grepl("^M_n_i$", colnames(dt)))==1){
    m_n <- dt$m_n_i
}

# Compute m_t
dt$m_t <- (dt$m_b * (p * dt$CG_CN  + (1-p) * CN_norm) - m_n  * (1-p) * CN_norm) / (p * dt$CG_CN)
# Note that %m_tumour = nA*mA + nB*mB for cases of allele specific methylation

# Correct methylation rate and calculate tumour-normal difference
dt$m_t_corr <- ifelse(dt$m_t<0,0,ifelse(dt$m_t>1,1,dt$m_t))
colnames(dt)[which(colnames(dt)=="m_t")] <- "m_t_raw"

# Caculate the coverage in the Tumour for each CpG cov * (p * CN / (p * CN + 2 * (1-p)))
dt[, cov_t := round2((cov_b * p * CG_CN) / (p * CG_CN + CN_norm * (1 - p)),digits = 0)]

# Remove sites with less than 3 reads from the tumour
dt <- dt[cov_t >= min_cov_t, ]

# Calculate counts of 1s in the tumour and the coverage
# dt[, M_t := round2(m_t * cov_t,digits = 0)]
# dt[, UM_t := cov_t - M_t]

# Remove NAs (should not be any)
dt <- dt[!is.na(m_t_raw), ]

return(dt)
}

# Calculate HDI by simulation
#
# Computes highest density interval from a sample of representative values,
#   estimated as shortest credible interval for a unimodal distribution
# Arguments:
#' @param M_b counts methylated in the tumour
#' @param UM_b counts unmethylated in the tumour
#' @param M_n counts methylated in the normal
#' @param UM_n counts unmethylated in the normal
#' @param p tumour purity
#' @param CN total tumour copy number
#' @param CN_n total normal copy number
#' @param credMass default is 0.99
#' credMass is a scalar between 0 and 1, indicating the mass within the
#' credible interval that is to be estimated.
#' @return Value: HDIlim is a vector containing the limits of the HDI
HDIofMCMC = function(M_b,UM_b,M_n,UM_n,p,CN,CN_n,credMass=0.99) {
  
  # sampleVec is a vector of representative values for the methylation rate 
  # probability distribution at a CpG
  sampleVec1 <- rbeta(n = 2000, shape1 = M_b+1, shape2 = UM_b+1)
  sampleVec2 <- rbeta(n = 2000, shape1 = M_n+1, shape2 = UM_n+1)
  
  # get the mb and mn constant
  cons1 <- p*CN + CN_n*(1-p)
  cons2 <- CN_n*(1-p)
  
  # obtain the mt sample vector
  sampleVec_scaled <-  (cons1*sampleVec1 - cons2*sampleVec2)/(p*CN)
  
  # sort the mt sample vector
  sortedPts <- sort(sampleVec_scaled)
  
  # get the width of the nth percentile where n=credMass
  ciIdxInc <- ceiling(credMass * length(sortedPts)) 
  
  # get the diff between all pairs at a suitable width
  ciWidth <- diff(x = sortedPts, lag = ciIdxInc) 
  
  # calculate mt HDI
  HDIlim <- c(lower = sortedPts[which.min(ciWidth)], 
              upper = sortedPts[which.min(ciWidth)+ciIdxInc])
  
  return(HDIlim)
}


# Plot 2d density comparing methylation rates
# 
# Arguments:
#' @param dt Data table with methylation information per CpG
#' @param path Character path variable pointing to the desired working directory.
#' This is where the output will be stored and should be constant for all CAMDAC functions.
#' @return NULL

plot_2d_density <- function(dt, path){

  # get relevant colummn
  tmp <- dt[, .SD, .SDcols=c("m_b","m_t_corr")]
  if("M_n" %in% colnames(dt)){
      tmp[, "m_n" := dt[,.(m_n)]]
  }
  if("M_n_i" %in% colnames(dt)){
     tmp[, "m_n" := dt[,.(m_n_o)]]
  }

  # plot 2d density
  p1 <- ggplot(tmp, aes(x=m_n, y=m_b)) +  theme_classic() + 
       geom_bin2d(bins=50) + ggtitle("Before CAMDAC") +
       scale_fill_gradient2(trans = "log10", low = "white", mid = "yellow" , high = "red3", 
                             midpoint = 1.8, breaks = c(10, 100, 1E3, 1E4, 1E5), 
                             labels = parse(text=c("10^1","10^2","10^3","10^4","10^5"))) + 
       geom_text(data=data.frame(x=0.2, y=0.95, tex=paste0("cor=",tmp[,round2(cor(m_n, m_b), digits=3)])), 
                  aes(x=x, y=y, label=tex)) +
       guides(fill=guide_legend(title=expression(count~(log[10])))) +
       theme(legend.position="none") +
       ylab("bulk methylation rate") + xlab(paste0("normal methylation rate"))

  p2 <- ggplot(tmp, aes(x=m_n, y=m_t_corr)) +  theme_classic() + 
       geom_bin2d(bins=50) + ggtitle("After CAMDAC") + 
       scale_fill_gradient2(trans = "log10", low = "white", mid = "yellow" , high = "red3", 
                             midpoint = 1.8, breaks = c(10, 100, 1E3, 1E4, 1E5), 
                             labels = parse(text=c("10^1","10^2","10^3","10^4","10^5"))) + 
       geom_text(data=data.frame(x=0.2, y=0.95, tex=paste0("cor=",tmp[,round2(cor(m_n, m_t_corr), digits=3)])), 
                  aes(x=x, y=y, label=tex)) +
       guides(fill=guide_legend(title=expression(count~(log[10])))) +
       ylab("pure tumour methylation rate") + xlab(paste0("normal methylation rate"))

  #p3 <- ggplot(tmp, aes(x=m_b, y=m_t_corr)) +  theme_classic() + 
  #     geom_bin2d(bins=50) + 
  #     scale_fill_gradient2(trans = "log10", low = "white", mid = "yellow" , high = "red3", 
  #                          midpoint = 1.8, breaks = c(10, 100, 1E3, 1E4, 1E5), 
  #                           labels = parse(text=c("10^1","10^2","10^3","10^4","10^5"))) + 
  #                aes(x=x, y=y, label=tex)) +
  #     geom_text(data=data.frame(x=0.2, y=0.95, tex=paste0("cor=",tmp[,round2(cor(m_b, m_t_corr), digits=3)])), 
  #     guides(fill=guide_legend(title=expression(count~(log[10])))) +
  #     ylab("pure tumour methylation rate") + xlab(paste0("bulk tumour methylation rate"))
 
  gglist1 <- list(p1)
  gglist2 <- list(p2)
  my_layout <- matrix(c(rep(1,2), rep(2, 3)), nrow=1)
  gd <- arrangeGrob(grobs = c(gglist1, gglist2), layout_matrix = my_layout)
  Sys.sleep(1)
  ggsave(gd, filename = paste0(path, "/tumour_versus_normal_methylomes.pdf"), device="pdf",
             width=6.5, height = 3, unit="in")
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
