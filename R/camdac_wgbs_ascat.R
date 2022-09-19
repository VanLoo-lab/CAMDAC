# ASCAT-WGBS

# Sort genomic loci
sort_genomic_dt <- function(dt, with_chr=F){
  if (with_chr){
    fact_levels=paste0("chr", c(1:22,"X","Y"))
  } else {
    fact_levels = c(1:22,"X","Y")
  }
  dt[, chrom:=factor(chrom, levels=fact_levels)]
  return(dt[order(chrom, POS)])
}


# Randomise BAF for ASCAT
randomise_BAF <- function(BAF){
  # Create a vector of randomly selected boolean values the length of input data
  selector = base::sample(c(TRUE,FALSE), length(BAF), replace=T, prob=c(0.5, 0.5))
  # Randomly assign or flip BAF values for a more uniform profile
  rBAF = data.table::fifelse(selector, BAF, 1-BAF)
  return(rBAF)
}


# Load a sample's SNP profile from allele_counts
load_snp_profile <- function(ac_file, loci_files){
  
  # Load data and subset to SNP sites and relevant columns
  snps <- data.table::fread(
      ac_file,
      select=c("CHR", "chrom", "POS", "total_counts", "total_depth", "ref", "alt", "BAF"))
  onames = names(snps) # save original columns for later
  
  # Filter to SNPs certified for allele counting using annotations in loci files
  ascat_pos = rbindlist(lapply(loci_files, function(x){
    load(x)
    df = loci_subset[(!is.na(loci_subset$CNA) & !is.na(loci_subset$ASCAT) & !is.na(loci_subset$POS))]
    seqlevelsStyle(df) = "NCBI"
    dt = data.table(chrom=as.character(seqnames(df)), POS=df$POS)
    return(dt)
  }))
  setkey(ascat_pos, chrom, POS)
  snps = ascat_pos[snps[!is.na(POS)], on=.(chrom, POS)] # Filters SNPs to positions in ASCAT_POS
    
  # Select SNP loci based on heading
  snps <- snps[!is.na(BAF),c("chrom", "POS", "total_counts", "total_depth", "ref", "alt", "BAF")]
  
  # Randomise BAF for downstream ASCAT analysis
  snps[, BAFr:=randomise_BAF(BAF)] 
  
  # Ensure chromosome field is properly formatted
  snps$chrom = as.character(snps$chrom)

  return(snps)
}

# Adds BAFr_N, LogR_n, total_depth_n
annotate_normal <- function(tsnps, nsnps, min_cov){
  # Ensure normal BAF randomised and LogR is 0
  nsnps[,BAFr:=randomise_BAF(BAF)]
  nsnps[,LogR:=0]
  
  # Add suffix to normal snps, subset columns and join to tumour snps data
  nsnps = nsnps[, .(chrom, POS, total_depth, total_counts, BAF, BAFr, LogR)]
  to_suffix = c("total_depth","BAF", "BAFr","LogR", "total_counts")
  setnames(nsnps, old=to_suffix, new=paste0(to_suffix, "_n"))
  setkey(nsnps, chrom, POS)
  
  tsnps = tsnps[nsnps, on=c("chrom","POS")]
  
  # TODO: Apply coverage filter from config or earlier in pipeline
  tsnps = tsnps[total_counts_n >= min_cov]
  
  return(tsnps)
}

# Calculate the tumour LogR from the tumour and normal sample
calculate_logr <- function(tsnps){
  
  #  Tumour LogR is the log2 normalised ratio of tumour and normal coverage
  #  at each position.
  calc_logr <- function(sample_cov, normal_cov){
  tumour_normal_coverage <- sample_cov / normal_cov
  mean_coverage_ratio <- mean(tumour_normal_coverage, na.rm=T)
  LogR <- log2(tumour_normal_coverage / mean_coverage_ratio)
  }
  
  # Calculate LogR and drop any NA sites. These sites have no coverage/SNP loci in the normal.
  tsnps[, LogR:=calc_logr(total_depth, total_depth_n)]
  tsnps <- tsnps[!is.na(LogR)]
  return(tsnps) 
}

# TODO: Time this function. Currently approx 8 minutes on a small sample
# TODO: Convert existing GC-content files to fst for faster read/write?
annotate_gc <- function(tsample, gc_refs, max_window=10000, n_cores=1){
  chrom=as.character(tsample$chrom)
  start=tsample$POS
  LogR=tsample$LogR
  
  # Fix chromosome to UCSC format
  chrom=data.table::fifelse(startsWith(chrom, "chr"), chrom, paste0("chr",chrom))
  
  # Create data table with LogR values
  dt = data.table::data.table(seqnames=chrom, start=start, end=start, LogR=LogR)
  data.table::setkey(dt, seqnames, start, end)
  
  # Get GC-LogR correlations for all window sizes below maximum
  logging::loginfo("Running GC correlation check")
  doParallel::registerDoParallel(cores=n_cores)
  gc_correlations = foreach::foreach(gc_file=gc_refs) %dopar% {
    # Get window size from filename
    window_size = as.numeric(gsub(".csv.gz","", regmatches(gc_file, regexpr("(\\d+).csv.gz",gc_file))))
    
    # Skip large windows. Here we correct for GC bias at the level of insert size
    if(window_size > max_window){return(list(gc_corr=NA, GC=NA, window=NA))}
    
    # Load data
    gcdf = data.table::fread(gc_file)[, .(seqnames, start, end, GC)]
    setkey(gcdf, seqnames, start, end)
    overlap=data.table::foverlaps(dt, gcdf)
    gc_corr = abs(cor(overlap$GC, overlap$LogR))
    
    return(list(gc_corr=gc_corr, GC=overlap$GC, window=window_size))
  }
  doParallel::stopImplicitCluster()
  

  logging::loginfo("GC correlation check complete")
  best_corr=gc_correlations[[which.max(sapply(gc_correlations, '[[', "gc_corr"))]]
  return(cbind(tsample, data.table(GC=best_corr$GC)))
}

annotate_repli <- function(tsample, repli_file){
  # Set variables for vectors
  chrom = tsample$chrom
  start = tsample$POS
  end = tsample$POS
  
  # Load replication timing data
  repli = data.table::fread(repli_file)
  names(repli)[names(repli)=="chromosome"] = "chrom"
  setkey(repli, chrom, start, end)
  
  # Create table for sample data, ensuring UCSC chromosome names
  dt = tsample[, .(chrom, start=POS, end=POS, LogR=LogR)]
  chrom=dt$chrom
  chrom=data.table::fifelse(startsWith(chrom, "chr"), chrom, paste0("chr",chrom))
  dt$chrom = chrom
  
  # Find the repliseq data nearest to each SNP
  repli_ranges = GRanges(seqnames=repli$chrom, ranges=IRanges(start=repli$start, end=repli$end))
  tumour_ranges = GRanges(seqnames=dt$chrom, IRanges(start=dt$start, end=dt$end))
  rep_match = repli[IRanges::nearest(tumour_ranges, repli_ranges)]
  names(rep_match) = sapply(names(rep_match), function(x) if(grepl("chrom|start|end",x)){gsub("$","_repli",x)}else{x})
  nearest_repli = cbind(rep_match, dt)
  
  # Get LogR correlation for best cell line
  cell_line_cols = !grepl("chrom|start|end|LogR", names(nearest_repli))
  correl = apply(nearest_repli[, ..cell_line_cols], MARGIN=2, FUN=function(x) abs(cor(x, nearest_repli$LogR)))
  best_line = names(which.max(correl))
  print(best_line)
  print(correl)
  
  # Combine with original dataframe and return
  result = cbind(tsample, data.table(repli=nearest_repli[[best_line]]))
  
  return(result)
}

spline_regress_logr <- function(LogR, GC, repli){
  # TODO: Why df=5?
  model = lm(LogR ~ splines::ns(x=GC, df=5, intercept=T) +
               splines::ns(x=repli, df=5, intercept=T),
             y= FALSE, model=FALSE, na.action = "na.exclude")
  return(model$residuals)
}

# Function from ASCAT/CAMDAC-RRBS
split_genome_WGBS = function(chrom, POS) {
  # TODO: Must be sorted and chrom must match POS. Easier to use a dataframe/table?
  
  # Convert chromosomes to numeric, including X and Y
  # suppressWarnings() used to stop warning that NAs introduced after coercion.
  # This is simply an effect of the way fcase handles the final condition. No NAs present.
  chrom = suppressWarnings(
    data.table::fcase(
    chrom == "X", 23,
    chrom == "Y", 24,
    # chrom that doesn't match condition is simply returned
    rep_len(TRUE, length(chrom)), as.double(chrom)
  ))
  
  # Identify large GAP regions and chromosome segment boundaries (1MB)
  # Diff goes pairwise through vector calculating differences. Which tells us where these are
  # and +1 required as it's actually the diff from the first element
  holesOver1Mb = which(diff(POS)>=1000000)+1 
  # Finds the indexes of chromosome borders
  chrBorders = which(diff(as.numeric(chrom))!=0)+1
  # Holes is a sorted list of indexes where the value preceding is a 1MB hole.
  holes = unique(sort(c(holesOver1Mb,chrBorders)))
  
  startseg = c(1,holes)
  endseg = c(holes-1,length(chrom))
  
  chr=lapply(seq(length(startseg)), function(x) startseg[x]:endseg[x])
  
  return(chr)
}

assign_genotypes <- function(BAF, as_logical=F){
  geno = data.table::fcase(
    BAF < 0.15, TRUE,
    BAF > 0.85, TRUE,
    default=FALSE
  )
  
  if(!as_logical){
    geno = factor(data.table::fcase(
      geno==TRUE, "Homozygous",
      geno==FALSE, "Heterozygous"
    ), levels=c("Homozygous", "Heterozygous"))
  }
  
  return(geno)
}

#' load_ascat_bc
#'
#' Create an ascat.bc object from input data vectors. Data must be sorted by genomic co-ordinate.
#'
#' @param logr_t Tumor LogR
#' @param baf_t Tumor BAF, ideally randomised for NGS.
#' TODO: Complete params
#'
#' @return List. An object containing the following fields:
#'     - Tumor_LogR
#'     - Tumor_BAF
#'     - Germline_LogR
#'     - Germline_BAF
#'     - SNPpos. dataframe of SNP positions with columns `chrom, POS`
#'     - chr. Genome segments, output from `split_genome`
#'     - samples. character vector of sample names
#'     - gender. The patient's sex. "XX" or "XY".
#'     - genotypes. The patient's Genotype profile
#'     - chrs. A vector of chromosome identifiers
#'     - ch. List. An element for each chrom with the SNPpos df indexes for SNPs belonging to that chromosome.
#'     
#' @keywords internal
#' @noRd
load_ascat_bc <- function(logr_t, baf_t, logr_n, baf_n, chrom, POS, samples, sex){
  
  # Build ch, a numeric for each chromosome
  chrom_names = c(1:22,"X","Y")
  ch = lapply(chrom_names, function(x) which(as.character(chrom)==x))
  names(ch)=chrom_names
  
  # Extract genotypes from normal BAF.
  genotypes = assign_genotypes(baf_n)
  
  # Split genome to create `chr`, a list of items, one per chromosome.
  # Function used depends on bsseq lib.
  chr = split_genome_WGBS(chrom, POS)
  
  ascat.bc = list(
    Tumor_LogR=data.frame(Tumor_LogR=logr_t),
    Tumor_BAF=data.frame(Tumor_BAF=baf_t),
    Germline_LogR=data.frame(Germline_LogR=logr_n),
    Germline_BAF=data.frame(Germline_BAF=baf_n),
    Tumor_LogR_segmented=NULL, Tumor_BAF_segmented=NULL,
    Tumor_counts=NULL, Germline_counts=NULL,
    SNPpos=data.frame(Chr=chrom, Position=POS), chr=chr,
    samples=paste(samples, sep="."), 
    chrs=c(1:22,"X","Y"), ch=ch, 
    gender=sex, sexchromosomes=c("X", "Y"),
    genotypes=data.frame(ggtypes=genotypes),
    failedarrays=NULL
  )
  
  names(ascat.bc$Tumor_LogR) = samples
  names(ascat.bc$Tumor_BAF) = samples
  names(ascat.bc$Germline_LogR) = samples
  names(ascat.bc$Germline_BAF) = samples
  
  return(ascat.bc)
}

#' @title bseq_bool
#' A helper function for filtering out reference and alternate
#' SNPs where bisulfite conversion cannot be distinguished.
#' @param ref 
#' @param alt 
#' @noRd
bseq_bool <- function(ref, alt){
  return(!(
    (ref == "C" & alt == "T") |
      (ref == "T" & alt == "C") |
      (ref == "A" & alt == "G") |
      (ref == "G" & alt == "A"))
  )
}

# rm_low_cov_singletons
# Remove low cov singletons. These are low-confidence SNPs that may cause ASCAT to produce small spurious segments.
# This function is ported from CAMDAC-RRBS
rm_low_cov_singletons = function(dt_sample_SNPs,min=3){
  
  # subselect relevant columns
  dt <- dt_sample_SNPs[,c("chrom","POS","total_counts")]
  
  # flag all SNPs with cov < 10
  low_cov <- dt$total_counts <= min
  
  # get neighbouring SNP index
  ranges <- -5:5
  idxs <- outer(which(low_cov==TRUE), ranges, `+`)
  # Creates a length(low_cov) x 11 matrix, where column #6 is the position of low-cov SNPs in our array,
  # while columns either side give the indexes of SNPs -5 and +6
  # NM: received error previously so filtering out negatives
  idxs <- pmax(idxs, 1); # pmax returns the maximum of the two values, so if 0 is higher for any item in the array, it's replaced
  
  # get neighbouring SNP genomic coordinates
  POSS <- suppressWarnings(matrix(dt$POS[idxs],ncol=length(ranges)))
  # Subtract each row from the position of the low coverage SNP
  POSS <-abs(POSS-dt[low_cov==TRUE]$POS)
  # Set low-cov SNPs greater than 1Mb away from neighbour to NA
  idxs[POSS>1E6] <- NA
  
  # Get the coverage at these loci (low cov snps and neighbours)
  covs <- matrix(dt$total_counts[idxs],ncol=length(ranges))
  # Get the average coverage across the region. 
  mean_covs <- rep(as.numeric(NA), length=nrow(dt))
  mean_covs[low_cov==TRUE] <- rowMeans(covs,na.rm=TRUE)
  # Note the vector of mean coverage is the length of original data
  
  # Flag and remove low coverage singletons on the following criteria:
  # 1) Low coverage AND there is an NA in the mean coverage values ## PROBLEM: Isn't this capturing chromosome boundaries?
  low_cov_na_mean = low_cov & is.na(mean_covs)
  # 2) Mean coverage is less than the minimum coverage
  mean_below_min = !is.na(mean_covs) & mean_covs < min
  low_cov_singleton = low_cov_na_mean | mean_below_min
  
  dt_sample_SNPs
  
  # Flag and remove low coverage singletons
  low_cov_singleton <- ifelse(low_cov==FALSE,FALSE,
                              ifelse(low_cov==TRUE&is.na(mean_covs),TRUE,
                                     ifelse(mean_covs>min,TRUE,FALSE)))
  dt_sample_SNPs <- dt_sample_SNPs[low_cov_singleton==FALSE,]
  
  # WARNING: NAs filtered during this function may include genome gaps as we are searching for
  # SNPs at a large distance from other SNPs.
  return(dt_sample_SNPs)
}


# Use BAF from normal allele counts file
use_external_normal_baf <- function(tumour, external_ac_file, config){
  
  stopifnot(fs::file_exists(external_ac_file))
  
  # Load TSNPs
  tsnps_file = CAMDAC::build_output_name(tumour, config, "tsnps")
  tsnps = data.table::fread(tsnps_file)
  tsnps$chrom = as.character(tsnps$chrom)
  
  # Create a backup of the original tSNPs file for record
  # This is instead of appending a column as normal SNPs overlap may remove loci
  fs::file_copy(tsnps_file, fs::path(tsnps_file, ext="initial"), overwrite=T)
  
  # Load external allele counts file
  ext_ac = data.table::fread(external_ac_file)
  stopifnot(
    all(c("#CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth") %in% names(ext_ac))
  )
  setnames(ext_ac, "#CHR", "chrom")
  ext_ac[, chrom:=gsub("chr","", chrom)]
  ext_ac = sort_genomic_dt(ext_ac)
  
  # Overlap loci and calculate BAF
  # Note: Good_depth used to simply select ext_ac rows that are not NA
  setkey(tsnps, chrom, POS)
  tsnps = tsnps[ext_ac, , on=.(chrom, POS)][!is.na(BAF) & !is.na(Good_depth)]
  tsnps[, BAF_n:=data.table::fcase(
    ref=="C" & alt=="T", Count_T/(Count_C+Count_T),
    ref=="C" & alt=="A", Count_A/(Count_C+Count_A),
    ref=="C" & alt=="G", Count_G/(Count_C+Count_G),
    ref=="G" & alt=="C", Count_C/(Count_G+Count_C),
    ref=="G" & alt=="T", Count_T/(Count_G+Count_T),
    ref=="G" & alt=="A", Count_A/(Count_G+Count_A),
    ref=="A" & alt=="T", Count_T/(Count_A+Count_T),
    ref=="A" & alt=="G", Count_G/(Count_A+Count_G),
    ref=="A" & alt=="C", Count_C/(Count_A+Count_C),
    ref=="T" & alt=="C", Count_C/(Count_T+Count_C),
    ref=="T" & alt=="A", Count_A/(Count_T+Count_A),
    ref=="T" & alt=="G", Count_G/(Count_T+Count_G)
  )]
  
  tsnps[, BAFr_n:=randomise_BAF(BAF_n)]
  
  tsnps[, `:=`(
    Count_A=NULL, Count_C=NULL, Count_G=NULL, Count_T=NULL, Good_depth=NULL
  )]
  
  tsnps = tsnps[total_depth_n >= 10]
  tsnps = tsnps[!is.na(BAF_n)]
  
  data.table::fwrite(tsnps, tsnps_file, compress="gzip")
  return(tsnps_file)
}

write_acf_and_ploidy_file <- function(tsnps, ascat.output, ascat.frag, sample_prefix, outdir){
  # Get genotypes for het/hom counts. TRUE for Hom and FALSE for Het
  genos = assign_genotypes(tsnps$BAFr, as_logical=T)
  
  fdata <- data.frame(
    ploidy=ascat.output$ploidy,
    ACF=ascat.output$aberrantcellfraction,
    num_het_SNPs_seg = nrow(ascat.frag$Tumor_LogR_segmented),
    num_hom_SNPs_seg = nrow(ascat.frag$Tumor_BAF_segmented[[1]]),
    num_het_SNPs_camdac = sum(!genos),
    num_hom_SNPs_camdac = sum(genos),
    median_depth_camdac = median(tsnps$total_depth),
    median_n_depth_camdac = median(tsnps$total_depth_n)
  )
  rownames(fdata)=NULL
  
  # Write to ouptut
  write.table(fdata, file=fs::path(outdir, paste0(sample_prefix, ".ACF.and.ploidy.txt")), sep="\t", row.names=F, col.names=T, quote=F)
}

# TODO: Refactor function to not use CAMDAC objects
run_ascat.m2 <-function(tumour, tsnps, penalty=200, outdir, rho_manual=NA, psi_manual=NA){

  # TODO: Make this an accessor function on the CAMDAC object so that it's consistent throughout code
  sample_prefix = paste(tumour$patient_id, tumour$sample_id, sep=".")
  
  # Load ASCAT object
  ascat.bc <- load_ascat_bc(
    logr_t = tsnps$LogR, baf_t = tsnps$BAFr,
    logr_n = tsnps$LogR_n, baf_n = tsnps$BAFr_n,
    chrom = tsnps$chrom, POS = tsnps$POS,
    samples = sample_prefix,
    sex = tumour$patient_sex
  )
  
  # Plot raw data
  # ascat.plotRawData(ascat.bc, img.dir=outdir, img.prefix=sample_prefix) # base ASCAT plotter
  ascat.m.plotRawData(ascat.bc, outdir=outdir)
  
  # Perform ASPCF segmentation
  gg = list(germlinegenotypes=matrix(assign_genotypes(ascat.bc$Germline_BAF, as_logical=T)))
  ascat.frag <- ASCAT::ascat.aspcf(ascat.bc, ascat.gg=gg, penalty=penalty,
                            out.dir=outdir, out.prefix=sample_prefix)

  # Plot segmented data
  # ascat.plotSegmentedData(ascat.frag, img.dir=outdir, img.prefix=sample_prefix) # base ASCAT plotter
  ascat.m.plotSegmentedData(ascat.frag, fname=sample_prefix, outdir=outdir)
  
  # Run ASCAT
  ascat.output <- ASCAT::ascat.runAscat(ascat.frag, gamma = 1, img.dir=outdir, img.prefix=sample_prefix,
                                        rho_manual=rho_manual, psi_manual=psi_manual)

  # Write ACF and ploidy text file. Format adapted from CAMDAC-RRBS
  write_acf_and_ploidy_file(tsnps, ascat.output, ascat.frag, sample_prefix, outdir)
  
  # Return ASCAT results
  return(list(
    ascat.bc=ascat.bc,
    ascat.output=ascat.output,
    ascat.frag=ascat.frag))
}

# Winsorize extreme values in tumor BAF
# Our rule is that we only remove SNPs that fall on 0/1 and are outliers from the median
winsorize <- function(BAF){
  # Calculate the running median
  # param: k determins how far the running median will be calculated
  medianFilter <- function(x,k){
    n <- length(x)
    filtWidth <- 2*k + 1
    
    #Make sure filtWidth does not exceed n
    if(filtWidth > n){
      if(n==0){
        filtWidth <- 1
      }else if(n%%2 == 0){
        #runmed requires filtWidth to be odd, ensure this:
        filtWidth <- n - 1
      }else{
        filtWidth <- n
      }
    }
    
    runMedian <- stats::runmed(x,k=filtWidth,endrule="median")
    
    return(runMedian)
  }
  
  # Set data points to 
  psi <- function(d,z){
    # d is raw_baf - running_median. Set 
    xwin <- d
    # z is raw value for tau*SD of the MAD of the running median
    # If the difference is greater than z in either + or - direction, winsorize it to z
    xwin[d < -z] <- -z
    xwin[d > z] <- z
    return(xwin)
  }
  
  #Perform MAD winsorization:
  # Tau is how many SDs away from median MAD we winsorize
  # K is how many probes the running median is calculated against
  madWins <- function(x,tau=2.5,k=40,digits=4){
    # Calculate the running median median
    xhat <- medianFilter(x,k)
    # Get the difference and SD
    d <- x-xhat
    SD <- stats::mad(d)
    # Set the winsorization threshold, i.e. X standard deviations of the MAD
    z <- tau*SD 
    # xwin is the factor by which we adjust each value of xhat
    # If d is within our z range, we simply add it back to get the value of x, 
    # otherwise we add z to winsorize
    xwin <- xhat + psi(d, z)
    
    # Detect outliers, i.e. SNPs where winsorization has been applied
    # This is done by simply copying the winsorizing rules
    
    # Detect outliers
    outliers <- (d < -z) | (d > z)
    return(list(ywin=xwin,sdev=SD,outliers=outliers))
  }
  
  return(madWins(BAF))
}

