##  Obtain copy number, purity and plot data
##  Latest update: 04/05/2021
##  Version 1.0.0
##  Author: Elizabeth Larose Cadieux

#' Obtain allele-specific copy number profiles, tumour purity and plot SNP data
#' 
#' \code{run_ASCAT.m}
#'
#' @param patient_id Character variable containting the patient id number
#' @param sample_id Character variable with the (control or tumour) sample_id
#' @param sex Character variable with  the patient expressed as "XX" for female
#' or "XY" for male.
#' This is important for copy number profiling. If sex is unknown, put "XY" for now,
#' then look at the allelic imbalance (BAF) on X in the germline outside pseudo-
#' autosomal regions. If there are little to no heterozygous SNPs, the sample is likely male.
#' @param patient_matched_normal_id Character variable with the sample ID of the matched normal control
#' @param path Character path variable pointing to the desired working directory.
#' This is where the output will be stored
#' IMPORTANT: The function output directory will be the in the path variable working
#' directory under "./Copy_number/sample_id/".
#' @param path_to_CAMDAC Character variable containting the path to the CAMDAC dir
#' including dir name (e.g. "/path/to/CAMDAC/").
#' @param build Character variable corresponding to the reference genome used for alignment.
#' CAMDAC is compatible with "hg19", "hg38", "GRCH37","GRCH38".
#' @param min_tumour Numerical value correspdonding to the minimum counts in the tumour
#' sample for germline SNPs to be included (default:10)
#' @param min_normal Numerical value correspdonding to the minimum counts for germline
#' SNPs to be included (default:1)
#' @param n_cores Numerical value correspdonding to the number of cores for parallel processing
#' @param reference_panel_coverage Path to the reference panel for the coverage.
#'
#' @return Three text files with all the CpG loci and their SNP and/or CpG methylation info 

run_ASCAT.m <- function (patient_id,sample_id,sex,
                         patient_matched_normal_id=NULL,
                         path,path_to_CAMDAC,build,
                         min_normal=10,min_tumour=1,
                         n_cores=1, reference_panel_coverage=NULL){
 
  if(getOption("scipen")==0){options(scipen = 999)}
  normal_id <- ifelse(is.null(patient_matched_normal_id), "none", patient_matched_normal_id)
  
  # Load SNP data
  path_patient <- file.path(path, patient_id)
  f_nm = paste0(path_patient,"/Allelecounts/", sample_id, "/",
                patient_id, ".", sample_id, ".SNPs.CpGs.all.sorted.RData")
  load(f_nm)
  rm(f_nm)
  
  # store relevant column names
  cols <- c("chrom", "POS", "total_counts", "total_depth", "ref", "alt", "BAF")
  
  # Rename R object (normal samples)
  if(sample_id==normal_id){
    if(is.null(reference_panel_coverage)){
    # Annotate SNP type
    dt_normal[, type := factor(ifelse(BAF < 0.15,"Homozygous",
                               ifelse(BAF > 0.85,"Homozygous",
                               ifelse(BAF >= 0.15 & BAF <= 0.85, "Heterozygous",
                               NA))), 
                             levels = c("Homozygous", "Heterozygous"))]
      cols <- c(cols, "type")
    }
    dt_sample <- dt_normal
    rm(dt_normal)
  }
  
  # Rename R object (tumour samples)
  if(sample_id != normal_id){
    dt_sample <- dt_tumour
    rm(dt_tumour)
  }
  
  # Create output directory
  path_output<-paste0(path_patient, "/Copy_number/",sample_id,"/")
  dir.create(path_output, recursive=TRUE, showWarnings = FALSE)
  orig_dir <- getwd()
  
  # Set reference human genome build variables
  cat(paste("Data with build ", build, sep = " "), "\n", sep = "")
  if(build=="GRCH37"){build="hg19"} # set build to to assembly version disregarging UCSC/Ensembl
  if(build=="GRCH38"){build="hg38"}
  
  # Fromat chromosome names
  chr_names=c(1:22,"X","Y")
  x <- substr(as.character(dt_sample$chrom[1]),1,3)
  if(x=="chr"){dt_sample[, chrom := substr(chrom, 4, 5)]}
  rm(x)

  # Select relevant SNP columns from the allele counts output
  dt_sample_SNPs <- dt_sample[, .SD, .SDcols=cols] 
  dt_sample_SNPs <- dt_sample_SNPs[!is.na(BAF)&total_counts>0&chrom!="Y",]

  # Remove duplicates and rename rows
  # Duplicates may arise from heterozygous overlapping CCGGG/CCCGG CCG/CGG  
  dt_sample_SNPs <- dt_sample_SNPs[!duplicated(paste(chrom, POS, sep="_")),]
  rownames(dt_sample_SNPs) <- c(1:nrow(dt_sample_SNPs))
  
  # Load germline tag
  if(sample_id != normal_id){
    if(is.null(reference_panel_coverage)){

      germline_snps_f_name <- paste0(path_patient,"/Copy_number/",normal_id,"/",patient_id, 
                            ".", normal_id, ".SNPs.RData")

      load(germline_snps_f_name); rm(germline_snps_f_name)
      dt_normal_SNPs <- dt[,.SD, .SDcols=c(cols,"type")]
      rm(dt)

      # format
      idxs <- match(c("total_counts","total_depth","BAF"), colnames(dt_normal_SNPs))
      colnames(dt_normal_SNPs)[idxs] <- c("total_counts_n","total_depth_n","BAF_n")
      rm(idxs)
    } else {
       dt_normal_SNPs <- read_fst(reference_panel_coverage, as.data.table=TRUE)
       dt_normal_SNPs <- 
           dt_normal_SNPs[, .SD, .SDcols= c("chrom","POS","ref","alt","total_depth_n")]
       dt_normal_SNPs[, BAF_n := as.numeric(NA)]
       dt_normal_SNPs <- dt_normal_SNPs[total_depth_n >= min_normal,]
    }
    # double check chromosome format
    x <- substr(as.character(dt_normal_SNPs$chrom[1]),1,3)
    if(x=="chr"){dt_normal_SNPs[, chrom := substr(chrom, 4, 5)]}
    rm(x)

    # Order seqlevels
    dt_normal_SNPs[, chrom := factor(chrom, levels=chr_names, ordered = TRUE)]
    dt_sample_SNPs[, chrom := factor(chrom, levels=chr_names, ordered = TRUE)]

    # Rename rows
    rownames(dt_normal_SNPs) <- 1:nrow(dt_normal_SNPs)
    rownames(dt_sample_SNPs) <- 1:nrow(dt_sample_SNPs)

    # Overlap normal and tumour
    dt_sample_SNPs <- merge(dt_sample_SNPs, dt_normal_SNPs, by = c("chrom","POS","ref", "alt"))
    rm(dt_normal_SNPs)
    cat("Germline SNP info loaded succesfully!\n")
  } ; rm(cols)
  
  # Obatin SNP genotype from bulk if there is no patient-matched normal
  if(!is.null(reference_panel_coverage) & sample_id != normal_id){
     dt_sample_SNPs[, tmp := ifelse(BAF<0.5, (1-BAF)*total_counts, BAF*total_counts)]
     dt_sample_SNPs[, tmp := round2(tmp, digits =0)]

     # probabilistic approach to assign heterozygous SNPs directly from tumour BAF profiles
     is_het <- function(x, y, pbin=0.01, probHom=.99, na.rm=TRUE) {
       flag <- logical(length=length(x))
       flag <- !(pbinom(unlist(x),size=unlist(y),prob=probHom,log.p=FALSE)>pbin)
       flag <- ifelse(flag==TRUE, "Heterozygous", "Homozygous")
       return(flag)
     }

     dt_sample_SNPs[, type := is_het(tmp, total_counts)]
     dt_sample_SNPs[, type := factor(type, levels = c("Homozygous", "Heterozygous"))]
     dt_sample_SNPs[, tmp := NULL]
  }
  
  # Apply coverage treshold
  if(sample_id == normal_id){
     dt_sample_SNPs <- dt_sample_SNPs[total_counts >= min_normal,]
     min <- min_normal
     cat("Minimum counts treshold in the matched normal",min_normal,"count(s)\n")
  }
  
  if(sample_id != normal_id){
     dt_sample_SNPs <- dt_sample_SNPs[total_counts >= min_tumour,]
     min <- min_tumour
     cat("Minimum counts treshold in tumour sample set to",min,
         "count(s) \nMinimum counts treshold in the matched normal",min_normal,"count(s)\n")
  }
  
  # Add SNP loci ids
  dt_sample_SNPs[, BAFloci := c(1:nrow(dt_sample_SNPs))]

  # Create vector for randomised BAF
  n = nrow(dt_sample_SNPs)
  if(sample_id != normal_id){dt_sample_SNPs[, rBAF_n := vector(length=n, mode="numeric")]}
  dt_sample_SNPs[, rBAF := vector(length=n, mode="numeric")]
  
  # randomise A and B alleles
  selector = round(runif(n))
  if(sample_id != normal_id & is.null(reference_panel_coverage)){
    dt_sample_SNPs[which(selector==0), "rBAF_n" := BAF_n]
    dt_sample_SNPs[which(selector==1), "rBAF_n" := 1 - BAF_n]
    # dt_sample_SNPs[, table(BAF_n==rBAF_n)]
  }
  if(sample_id != normal_id & !is.null(reference_panel_coverage)) {
    dt_sample_SNPs[, rBAF_n := as.numeric(NA)]
  }
  dt_sample_SNPs[which(selector==0), "rBAF" := BAF]
  dt_sample_SNPs[which(selector==1), "rBAF" := 1 - BAF]
  rm(selector,n)
  
  # Whilst germline heterozygous SNPs can have LOH and appear homozygous in the tumour, 
  # Germline homozygous loci should remain homozygous in the tumour (par from mutations e.g. at CCT)
  # Calculate how many germline homozygous SNPs become het at various treshold to refine
  # heterozygous BAF treshold in RRBS
  
  if(sample_id!=normal_id){
    # Remove low coverage singletons
    cat("Removing low coverage singletons\n")
    n <- nrow(dt_sample_SNPs)
    dt_sample_SNPs <- remove_low_cov_singletons(dt_sample_SNPs=dt_sample_SNPs,min=min)
    cat(paste0("Low coverage singletons removed (", 
        round2((1-(nrow(dt_sample_SNPs)/n))*100, digits=2),"% of SNPs).\n"))
    rm(n)

    # Set reference file names for LogR bias correction
    if(is.null(reference_panel_coverage)){
      fragments_file=paste0(path_patient,"/Copy_number/",normal_id,"/msp1_fragments_RRBS.RData")
    } else {
      fragments_file=paste0(path_patient,"/Copy_number/",sample_id,"/msp1_fragments_RRBS.RData")
    }
    replic_timing_file_prefix <- 
    paste0(file.path(path_to_CAMDAC, "pipeline_files"),
           "/1000genomes_2012_v3_repliTiming/1000_genomes_replication_timing_chr_")
  
    # Compute and correct logR
    dt_stats <- LogR_correction(dt_sample=dt_sample,dt_SNPs=dt_sample_SNPs,
                                build=build,chr_names=chr_names,
                                min_normal=min_normal,
                                fragments_file=fragments_file,
                                replic_timing_file_prefix=replic_timing_file_prefix,
                                n_cores=n_cores)
    rm(fragments_file,replic_timing_file_prefix)  
    # Note that positions w/ undetermined fragment size or GC content will be removed
  
    # format output for overlap with sample SNP data.table
    dt_stats[, POS := start] 
    dt_stats[, c("start", "end") := list(NULL, NULL)]
    setkeyv(dt_sample_SNPs, c("chrom","POS"))
    setkeyv(dt_stats, c("chrom","POS"))
    dt_sample_SNPs <- dt_sample_SNPs[dt_stats, nomatch=0]
    rm(dt_stats)
  
    cat("LogR correction completed\n")
  } else {      
    # format normal seqnames in normal
    y <- substr(as.character(dt_sample_SNPs$chrom[1]),1,3)
    if(y=="chr"){dt_sample_SNPs[, chrom := substr(chrom, 4, 5)]}
    rm(y)
  }
  
  # order the dt
  dt_sample_SNPs[, chrom := factor(chrom, levels=chr_names, ordered = TRUE)] 
  dt_sample_SNPs <- dt_sample_SNPs[order(chrom, POS)]
  dt_sample_SNPs[, BAFloci := 1:nrow(dt_sample_SNPs)]
  
  # Select columns
  cols <- c("chrom","POS","total_counts","total_depth","ref","alt","BAF","rBAF")
  if(sample_id!=normal_id){
    cols=c(cols,"LogR_t","LogR_t_corr","total_depth_n","total_counts_n",
           "BAF_n","rBAF_n","LogR_n","msp1_length","GC_content","replic")
    if(!is.null(reference_panel_coverage)){cols <- cols[-which(cols=="total_counts_n")] }
  }
  dt <- dt_sample_SNPs[, .SD, .SDcols=c(cols,"type")]
 
  # save dt
  save(dt, file=paste(path_output, patient_id, ".", sample_id, ".SNPs.RData", sep = ""))
  rm(dt)

  # Switch WD for ASCAT
  setwd(path_output)

  if(sample_id != normal_id){
    # Annotate
    dt_sample_SNPs[, chrom := as.character(chrom)]
    SNPpos = dt_sample_SNPs[, .SD, .SDcols=c("chrom","POS","BAFloci")]
    ch = list()
    for (i in 1:length(chr_names)) {
      temp = which(as.character(SNPpos$chrom)==chr_names[i])
      if (length(temp) == 0) {
        ch[[i]] = 0
      } else {
        ch[[i]] = temp[1]:temp[length(temp)]
      }
    }
    
    chr=split_genome_RRBS(SNPpos)
    types <- dt_sample_SNPs[, as.character(type)=="Homozygous"]
    
    ascat.bc = list(Tumor_LogR=data.frame(LogR_t_corr=dt_sample_SNPs$LogR_t_corr),
                    Tumor_BAF=data.frame(rBAF_t=dt_sample_SNPs$rBAF),
                    Germline_LogR=data.frame(LogR_n=dt_sample_SNPs$LogR_n),
                    Germline_BAF=data.frame(rBAF_n=dt_sample_SNPs$rBAF_n),
                    Tumor_LogR_segmented=NULL, Tumor_BAF_segmented=NULL,
                    Tumor_counts=NULL, Germline_counts=NULL,
                    SNPpos=SNPpos[,c("chrom","POS")], chr=chr,
                    samples=paste(patient_id, sample_id, sep="."), 
                    chrs=chr_names, ch=ch, 
                    gender=sex, sexchromosomes=c("X", "Y"),
                    genotypes=data.frame(ggtypes=types))
    rm(chr,ch,SNPpos)

    # chrs = chromosome_names
    # ch = list of length chr_names. Each list contrains all the position for any one of chr_names
    # chr = All SNP loci from 1:n split in different lists where gaps of more than 1Mb are
    # found or chromosome borders
    
    ascat.m.plotRawData(ascat.bc, raw_LogR=dt_sample_SNPs$LogR_t, pch = 10, cex = 0.2, lim_logR = 2.5)
    save(ascat.bc, file = paste(patient_id, sample_id, "ascat.bc.RData", sep = "."))
    cat("ASCAT object created\n")
    
    # Carry out segmentation
    gg = list(germlinegenotypes=ascat.bc$genotypes)
    ascat.frag <- ascat.aspcf(ascat.bc, ascat.gg=gg, penalty=200)
    # penalty = 200 recommended for sequencing data

    # fix issue with ascat.ascpcf renaming samples 
    ascat.frag$samples <- paste(patient_id, sample_id, sep=".")
    rm(ascat.bc)
    
    ascat.m.plotSegmentedData(ascat.frag, lim_logR = 2.5) 
    save(ascat.frag, file = paste(patient_id, sample_id, "ascat.frag.RData", sep = "."))
    cat("\nASCAT copy number segmentation completed\n")
    
    # Run copy number caller a first time to get the distance matrix
    ascat.output <- ascat.runAscat(ascat.frag, gamma = 1)
    save(ascat.output, file = paste(patient_id, sample_id,"ascat.output.RData", sep = "."))
    num_het_SNPs = nrow(ascat.frag$Tumor_LogR_segmented)
    num_hom_SNPs = nrow(ascat.frag$Tumor_BAF_segmented[[1]])
    rm(ascat.frag)
    
    if(file.exists(paste(patient_id, sample_id,"ASCATprofile.png", sep = "."))){
      cat("\nASCAT completed\n")
    
      # Save purity and ploidy
      f.nm <- paste(patient_id, ".", sample_id,".ACF.and.ploidy.txt", sep = "")
      file.create(f.nm)
      f <- file(f.nm, open="w") 
      
      # Save/write
      dt <- data.frame(ploidy=ascat.output$ploidy,ACF=ascat.output$aberrantcellfraction, 
                      num_het_SNPs = num_het_SNPs, num_hom_SNPs = num_hom_SNPs,
                      #mean_depth = mean(dt_sample_SNPs$total_depth, na.rm=TRUE), 
                      median_depth = median(dt_sample_SNPs$total_depth, na.rm=TRUE),
                      median_n_depth = median(dt_sample_SNPs$total_depth_n, na.rm=TRUE))
      rm(ascat.output, num_het_SNPs, num_hom_SNPs)
      cat(format_delim(dt, delim = "\t", col_names = T),  file = f)
      close(f); rm(dt,f,f.nm)
      cat(paste("\nPloidy, Purity and summary stats saved in ",
                path_output,patient_id,".",sample_id,".ACF.and.ploidy.txt","\n",sep = ""))
      } else {
        cat("\nASCAT could not find a solution\n")  
    }
    
    # convert to data.table
    dt_sample_SNPs[, chrom := as.character(chrom)]

    # format seqnames
    x <- dt_sample[1, substr(as.character(CHR),1,3)]
    if(x=="chr"){dt_sample[, chrom := substr(dt_sample$CHR, 4, 5)]}
    if(x!="chr"){dt_sample[, chrom := dt_sample$CHR]}
    
    # join SNP and methylation info for plot
    dt_sample <- dt_sample[!is.na(POS),]
    dt_sample$chrom <- as.character(dt_sample$chrom)
    dt <- merge(dt_sample_SNPs, 
                dt_sample[,c("chrom","POS","ref","alt","BAF","m","width","CCGG")], 
                by=c("chrom","POS","ref","alt","BAF"))
    dt[, flag := factor(ifelse(CCGG>0,"CCGG",
                        ifelse(width>=2&!is.na(dt$m), "CG", 
                        "neither")),
                levels=c("CCGG", "CG", "neither"))]
    rm(dt_sample)
    
    # run plot function
    outfile = paste(patient_id, sample_id,"SNP_data.pdf", sep = "_")
    plot_SNP_info(dt=dt,outfile=outfile,min=min)
    cat("BAF and LogR diagnostics plots generated\n")
    }
    
    if(is.null(reference_panel_coverage)&sample_id==normal_id){ 
      # convert to data.table
      dt_sample_SNPs[, chrom := as.character(chrom)]

      # join SNP and methylation info for plot
      dt_sample <- dt_sample[!is.na(POS),]
      dt_sample$chrom <- as.character(dt_sample$chrom)
      dt <- merge(dt_sample_SNPs, 
                dt_sample[,c("chrom","POS","ref","alt","BAF","m","width","CCGG")], 
                by=c("chrom","POS","ref","alt","BAF"))
      dt[, flag := factor(ifelse(CCGG>0,"CCGG",
                        ifelse(width>=2&!is.na(dt$m), "CG", 
                        "neither")),
                levels=c("CCGG", "CG", "neither"))]
      rm(dt_sample)
    
      # set vars
      min=min_normal
      outfile=paste(patient_id, sample_id, "normal_SNP_data.pdf", sep = "_")
      
      # run plot function
      plot_normal_SNP_info(dt=dt,outfile=outfile,min=min)
      cat("Normal BAF plots generated\n")
    }

  setwd(orig_dir)

}

#' @title split_genome_RRBS
#' @author Peter Van Loo, modified by Elizabeth Larose Cadieux
#' @description A helper function to split the genome into parts
#' @param SNPpos A data.frame with a row for each SNP. 
#' The first column is chromosome, second column position
#' @noRd

split_genome_RRBS = function(SNPpos) {
  # look for gaps of more than 1Mb and chromosome borders
  SNPposnum <- SNPpos
  SNPposnum[, chrom := ifelse(chrom == "X", 23, chrom)]
  SNPposnum <- SNPposnum[, lapply(.SD, as.numeric)]
  
  # identify large regions chromosome and segment boundaries
  holesOver1Mb = which(diff(SNPposnum$POS)>=1000000)+1
  chrBorders = which(diff(as.numeric(SNPposnum$chrom))!=0)+1
  holes = unique(sort(c(holesOver1Mb,chrBorders)))
  
  # find which segments are too small
  joincandidates=which(diff(c(0,holes,dim(SNPposnum)[1]))<200)
  
  # if it's the first or last segment, just join to the one next to it, irrespective of chromosome and positions
  while (1 %in% joincandidates) {
    holes=holes[-1]
    joincandidates=which(diff(c(0,holes,dim(SNPposnum)[1]))<200)
  }
  while ((length(holes)+1) %in% joincandidates) {
    holes=holes[-length(holes)]
    joincandidates=which(diff(c(0,holes,dim(SNPposnum)[1]))<200)
  }
  
  while(length(joincandidates)!=0) {
    # the while loop is because after joining, segments may still be too small..
    startseg = c(1,holes)
    endseg = c(holes-1,dim(SNPposnum)[1])
    
    # for each segment that is too short, see if it has the same chromosome as the segments before and after
    # the next always works because neither the first or the last segment is in joincandidates now
    previoussamechr = SNPposnum[endseg[joincandidates-1],1]==SNPposnum[startseg[joincandidates],1] 
    nextsamechr = SNPposnum[endseg[joincandidates],1]==SNPposnum[startseg[joincandidates+1],1]
    
    distanceprevious = SNPposnum[startseg[joincandidates],2]-SNPposnum[endseg[joincandidates-1],2]
    distancenext = SNPposnum[startseg[joincandidates+1],2]-SNPposnum[endseg[joincandidates],2]
    
    # if both the same, decide based on distance, otherwise if one the same, take the other, if none, just take one.
    joins = ifelse(previoussamechr&nextsamechr, 
                   ifelse(distanceprevious>distancenext, joincandidates, joincandidates-1),
                   ifelse(nextsamechr, joincandidates, joincandidates-1))
    
    holes=holes[-joins]
    
    joincandidates=which(diff(c(0,holes,dim(SNPposnum)[1]))<200)
  }

  startseg = c(1,holes)
  endseg = c(holes-1,dim(SNPposnum)[1])
  
  chr=list()
  for (i in 1:length(startseg)) {
    chr[[i]]=startseg[i]:endseg[i]
  }
  
  return(chr)
}

#' @title remove_low_cov_singletons
#' @description Remove low coverage singletons outliers
#' @author Elizabeth larose cadieux
remove_low_cov_singletons = function(dt_sample_SNPs,min){

  # subselect relevant columns
  dt <- dt_sample_SNPs[,c("chrom","POS","total_counts")]
  
  # flag all SNPs with cov < 10
  low_cov <- dt$total_counts <= min
  
  # get neighbouring SNP index
  ranges <- -5:5
  idxs <- outer(which(low_cov==TRUE), ranges, `+`)
  idxs[idxs<0] <- 0
  
  # get neighbouring SNP genomic coordinates
  POSS <- matrix(dt$POS[idxs],ncol=length(ranges))
  POSS <- abs(POSS-dt[low_cov==TRUE, "POS"])
  idxs[POSS>1E6] <- NA
  
  # get neighbouring SNP genomic coordinates
  covs <- matrix(dt$total_counts[idxs],ncol=length(ranges))
  mean_covs <- rep(as.numeric(NA), length=nrow(dt))
  mean_covs[low_cov==TRUE] <- rowMeans(covs,na.rm=TRUE)

  # Flag and remove low coverage singletons
  low_cov_singleton <- ifelse(low_cov==FALSE,FALSE,
                       ifelse(low_cov==TRUE&is.na(mean_covs),TRUE,
                       ifelse(mean_covs>min,TRUE,FALSE)))
  dt_sample_SNPs <- dt_sample_SNPs[low_cov_singleton==FALSE,]
  
  return(dt_sample_SNPs)
  #rm(POSS,idxs,mean_covs,covs,ranges,dt,low_cov,low_cov_singleton)
}

#' @title LogR_correction
#' @description Correct logR for msp1 fragment size bias and GC content
#' @author Elizabeth Larose Cadieux 
#' @param dt_sample Allelecounts output as a data.table
#' @param dt_SNPs Allelecounts output subset to QC'ed SNP positions
#' @param build Character variable corresponding to the reference genome version used for alignment
#' @param chr_names Character variable with the seqlevels.
#' @param min_normal Numerical with the minimum normal coverage threshold
#' @param fragments_file CAMDAC reference MspI fragments file
#' @param replic_timing_file_prefix CAMDAC reference replication timing files path and file name prefix
#' @param n_cores Numerical value correspdonding to the number of cores for parallel processing

LogR_correction = function(dt_sample,dt_SNPs,build,chr_names,min_normal,
                           fragments_file,replic_timing_file_prefix,n_cores){
  
  # Load MspI fragments files
  load(fragments_file)
  
  # Format fragments GRanges as data.table
  fragments_patient <- data.table(fragments_patient)
  setkeyv(fragments_patient, c("CHR", "start", "end"))
  
  # Format SNP genomic coordinates
  dt_SNPs[, c("start", "end") := list(POS, POS)]
  x <- substr(as.character(dt_SNPs$chrom[1]),1,3)
  if(x!="chr"){dt_SNPs[, CHR := paste0("chr",chrom)]}
  if(x=="chr"){dt_SNPs[, CHR := chrom]}
  rm(x)
  
  # Format as data.table
  dt_SNPs <- dt_SNPs[,c("CHR","start","end","total_depth","total_depth_n")]
  setkeyv(dt_SNPs, c("CHR", "start", "end"))

  # We only need the fragments overlapping with SNP positions
  overlaps <- foverlaps(fragments_patient, dt_SNPs, nomatch=0)
  
  # Save SNP with fragment ID and depth for later LogR correction
  dt_SNPs <- overlaps[!is.na(total_depth_n), .SD,
                      .SDcols=c("CHR","start","end","total_depth","total_depth_n",
                                "msp1_length","ID")]
  # Duplicates created due to SNPs overlaps with 2 fragments (i.e. SNP at a CCGG)
  # Keep for now, will take mean logR
  
  # Keep fragment that overlap with SNPs
  fragments_patient <- overlaps[!is.na(total_depth_n), .SD,
                                .SDcols=c("CHR","i.start","i.end","msp1_length",
                                          "CG","C","Th","A","G","ID")]
  fragments_patient <- fragments_patient[!duplicated(fragments_patient)]
  colnames(fragments_patient) <- gsub("i\\.", "", colnames(fragments_patient))
  setkeyv(fragments_patient, c("CHR", "start", "end"))
  rm(overlaps)
  
  # Get mean methylation rate for each fragment for GC content calculation
  dt_sample_m <- dt_sample[!is.na(m), .SD, .SDcols=c("CHR","start","end","m")];# Load CpGs
  dt_sample_m[,c("start","end","m")] <- lapply(dt_sample_m[,c("start","end","m")], as.numeric)
  
  # Get data per unique loci
  dt_sample_m <- dt_sample_m[,.(m = mean(m, na.rm = T)), keyby = .(CHR, start, end)]
  
  # Change coords of CCGG to CG
  dt_sample_m[, width := end-start+1]
  dt_sample_m[, start := ifelse(width == 4, start + 1, start)]
  dt_sample_m[, end := ifelse(width == 4, end - 1, end)]
  
  # Format data.table
  dt_sample_m <- dt_sample_m[, .SD, .SDcols=c("CHR","start","end","m")]
  setkeyv(dt_sample_m, c("CHR", "start", "end"))
  
  # Overlap methylation data with fragments
  overlaps <- foverlaps(fragments_patient, dt_sample_m)
  dt <- overlaps[!is.na(m), .SD, .SDcols=c("CHR","i.start","i.end","msp1_length",
                                         "CG","C","Th","A","G","ID", "m")]
  colnames(dt) <- gsub("i\\.", "", colnames(dt))
  rm(overlaps)
  
  # Get the methylation rate per fragment
  dt_SNPs_m <- dt[,.(mean_m = mean(m, na.rm = T)),
                      keyby = .(CHR, start, end, msp1_length, ID, CG, C, Th, A, G)]
  rm(dt)
  
  # Add fragments with SNP but no methylation info
  IDs <- !fragments_patient$ID%in%dt_SNPs_m$ID
  dt_SNPs_m2 <- as.data.frame(fragments_patient[IDs]);rm(IDs)
  
  mean_m1 <- ifelse((dt_SNPs_m2$ID-2)%in%dt_SNPs_m$ID, dt_SNPs_m$mean_m[match(dt_SNPs_m2$ID-2,dt_SNPs_m$ID)],NA)
  mean_m2 <- ifelse((dt_SNPs_m2$ID-1)%in%dt_SNPs_m$ID, dt_SNPs_m$mean_m[match(dt_SNPs_m2$ID-1,dt_SNPs_m$ID)],NA)
  mean_m3 <- ifelse((dt_SNPs_m2$ID+1)%in%dt_SNPs_m$ID, dt_SNPs_m$mean_m[match(dt_SNPs_m2$ID+1,dt_SNPs_m$ID)],NA)
  mean_m4 <- ifelse((dt_SNPs_m2$ID+2)%in%dt_SNPs_m$ID, dt_SNPs_m$mean_m[match(dt_SNPs_m2$ID+2,dt_SNPs_m$ID)],NA)
  mean <- mapply(function(mean_m1,mean_m2,mean_m3,mean_m4){
                  mean<-ifelse(sum(!is.na(cbind(mean_m1,mean_m2,mean_m3,mean_m4)))>0,
                  mean(cbind(mean_m1,mean_m2,mean_m3,mean_m4),na.rm=T),0.5)# Assume m = 0.5 when no other info 
                 },mean_m1,mean_m2,mean_m3,mean_m4)

  dt_SNPs_m2$mean_m <- mean
  rm(list=ls(pattern="mean"))
  
  # Join data.tables with methylation info
  dt_SNPs_m <- rbind(dt_SNPs_m, dt_SNPs_m2);rm(dt_SNPs_m2)
  dt_SNPs_m <- dt_SNPs_m[order(factor(CHR, levels=paste0("chr", c(1:22,"X","Y")), ordered=TRUE), start, end)]
  
  # Calculate GC content, correcting for CG methylation
  dt_SNPs_m[, C_read := mean_m*CG]
  # could use all_counts to weigh each contributing meth rates
  dt_SNPs_m[, T_read := dt_SNPs_m[,(C - mean_m*CG + Th)]]
  
  # Get GC content (mean human genome GC content is ~ 46%)
  dt_SNPs_m[, GC_content := dt_SNPs_m[,(C_read+G)/(C_read+G+A+T_read)]]
  
  # Annotate GC for each SNP 
  index <- dt_SNPs_m[,(!duplicated(ID,GC_content))]
  dt_SNPs <- dt_SNPs[dt_SNPs_m[index,.(ID,GC_content)], on="ID", nomatch=0] 
  rm(index, dt_SNPs_m)
  
  # Get data per unique loci
  assign("tmp",dt_SNPs)
  dt_SNPs <- tmp[,.(msp1_length=mean(msp1_length, na.rm = T), 
                        GC_content=mean(GC_content, na.rm = T),
                    ID=ID[which.max(msp1_length)]),
                     keyby=.(CHR,start,end,total_depth,total_depth_n)]
  rm(tmp)
  
  # Get frag_size
  dt_SNPs[, Log_msp1_length := log10(msp1_length)]
  
  # Calculate logR
  # assume that normallogR is 0, and normalise mutantLogR to normalLogR
  dt_SNPs$LogR_n <- vector(length=nrow(dt_SNPs), mode="integer") 
  dt_SNPs$LogR_t <- dt_SNPs[,(total_depth / total_depth_n)]
  tmp <- dt_SNPs[,log2(LogR_t/mean(LogR_t, na.rm=TRUE))]
  dt_SNPs$LogR_t <- tmp
  rm(tmp)
  
  # load replication timing files 
  chrom_idx = 1:23
  if(build=="hg19"){replic_files = paste0(replic_timing_file_prefix, chrom_idx, ".fst")}
  if(build=="hg38"){replic_files = paste0(replic_timing_file_prefix, chrom_idx,"_",build,".fst")}
  replic_data = data.table(do.call(rbind, mclapply(replic_files, fst::read_fst, mc.cores = n_cores)))
  
  # format replication timing data
  cols <- colnames(replic_data)
  replic_data[, 2:ncol(replic_data)] <- replic_data[, lapply(.SD, as.numeric), .SDcols= 2:ncol(replic_data)]
  replic_data[, chr := factor(chr, levels=c(1:22, "X", "Y"), ordered = TRUE)]
  replic_data <- replic_data[order(chr, pos),]
  rm(replic_files)
  
  # check seqlevels match
  colnames(replic_data)[1] <- "CHR"
  x <- substr(replic_data$CHR[1], 1, 3)=="chr"
  if(x==FALSE){
      replic_data[, CHR := paste0("chr", CHR)]
  } 
  replic_data[, c("start", "end") := list(pos, pos)]
  replic_data[, c("pos") := NULL ]
  rm(x)
  
  # extend to 1 kb window either side of replic data points
  replic_data[, start := start - 1000]
  replic_data[, end := end + 1000]

  # Annotate fragments with replic data
  colnames(fragments_patient) <- gsub("i.", "", colnames(fragments_patient))
  setkeyv(replic_data, c("CHR", "start", "end"))
  overlaps <- foverlaps(replic_data, fragments_patient, nomatch = NA)
  overlaps <- overlaps[!is.na(ID)]
   
  # annotate each fragment with replic_data
  # take the mean replic timing value per fragment
  cols_ids <- (which(colnames(overlaps)=="ID")+1):(which(colnames(overlaps)=="i.start")-1)
  fragments_replic_data <- overlaps[, lapply(.SD, mean, na.rm=TRUE), .SDcols=colnames(overlaps)[cols_ids], by=.(ID)]
  rm(overlaps, replic_data, cols_ids)

  locimatches <- match(x = dt_SNPs[, ID], table = fragments_replic_data[,ID])
  dt_SNPs <- dt_SNPs[!is.na(locimatches), ]
  fragments_replic_data <- fragments_replic_data[na.omit(locimatches), ]
  
  # Find and select cell line Repli-Seq dataset with highest correlation to the observed LogR biases 
  corr_rep = abs(cor(fragments_replic_data[, .SD, .SDcols=2:ncol(fragments_replic_data)], 
                     dt_SNPs$LogR_t, use="complete.obs")[,1])
  maxreplic = which.max(corr_rep)+1 ;rm(corr_rep)
  cat(paste0("Replication timimg correction based on ", colnames(fragments_replic_data)[maxreplic], 
             " ENCODE cell line Repli-Seq data."))
  
  # annotate each SNP with fragment replication timing info
  dt_SNPs$replic <- fragments_replic_data[, .SD, .SDcols=maxreplic]
  rm(maxreplic,fragments_replic_data,locimatches)
  
  # Remove low cov sites in the normal 
  # (low cov in the normal could be losses, low cov in normal are problematic to sequence, not abberations)
  dt_SNPs <- dt_SNPs[total_depth_n >= min_normal,] # should already be removed in previous steps
  
  # Combining GC and frag size analysis in a linear model
  dt_SNPs <- dt_SNPs[order(CHR, start, end),]
  dt_SNPs <- dt_SNPs[!duplicated(dt_SNPs),]
  
  y<-dt_SNPs$LogR_t
  x1<-dt_SNPs$Log_msp1_length
  x2<-dt_SNPs$GC_content
  x3<-dt_SNPs$replic
  
  model = lm(y ~ splines::ns(x = x1, df = 5, intercept = TRUE) + 
                 splines::ns(x = x2, df = 5, intercept = TRUE) + 
                 splines::ns(x = x3, df = 5, intercept = TRUE), y=FALSE, 
                 model = FALSE, na.action="na.exclude")
  #model<-lm(y~poly(x1,2)+poly(x2,2),y=TRUE)#;plot(model, pch=19,cex=0.02)
  dt_SNPs[, LogR_t_corr := model$residuals] ; rm(model)

  # format seqnames
  colnames(dt_SNPs)[1] <- "chrom"
  x <- substr(dt_SNPs$chrom[1], 1, 3)=="chr"
  if(x==TRUE){
    dt_SNPs[, chrom := substr(chrom, 4, 5)]
  }   

  # return data.table with corrected LogR and coveriates
  return(dt_SNPs)
}

#' @title ascat.m.plotSegmentedData
#' @description Plot segmentated BAF LogR 
#' @param ASCATobj an ASCAT object (e.g. data structure from ascat.loadData)
#' 
#' @return Produces png files showing the logR and BAF values for tumour and germline samples
#' @author Peter Van Loo
#' @export

ascat.m.plotSegmentedData <- function (ASCATobj, lim_logR=2) 
{
  for (arraynr in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    Select_nonNAs = rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]])
    AllIDs = 1:dim(ASCATobj$Tumor_LogR)[1]
    names(AllIDs) = rownames(ASCATobj$Tumor_LogR)
    HetIDs = AllIDs[Select_nonNAs]
    png(filename = paste(ASCATobj$samples[arraynr], ".ASPCF.png", 
                         sep = ""), width = 2000, height = 1000, res = 200)
    par(mar = c(0.5, 5, 5, 0.5), mfrow = c(2, 1), cex = 0.4, 
        cex.main = 3, cex.axis = 2)
    r = ASCATobj$Tumor_LogR_segmented[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), 
                                      arraynr]
    beta = ASCATobj$Tumor_BAF_segmented[[arraynr]][, , drop = FALSE]
    plot(c(1, length(r)), c(-lim_logR ,lim_logR), type = "n", xaxt = "n", 
         main = paste(ASCATobj$samples[arraynr], ", LogR", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_LogR[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), 
                               arraynr], col = "red", pch = 10, cex = 0.20)
    points(r, col = "blue")
    abline(v = 0.5, lty = 1, col = "lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = intersect(ASCATobj$ch[[j]], HetIDs)
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len
      tpos = (chrk_tot_len + chrk_tot_len_prev)/2
      text(tpos, lim_logR-0.5, ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v = vpos + 0.5, lty = 1, col = "lightgrey")
    }
    plot(c(1, length(beta)), c(0, 1), type = "n", xaxt = "n", 
         main = paste(ASCATobj$samples[arraynr],", BAF", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_BAF[rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]), 
                              arraynr], col = "red", pch = 10, cex = 0.20)
    points(beta, col = "blue")
    points(1 - beta, col = "blue")
    abline(v = 0.5, lty = 1, col = "lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = intersect(ASCATobj$ch[[j]], HetIDs)
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len
      tpos = (chrk_tot_len + chrk_tot_len_prev)/2
      text(tpos, 1, ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v = vpos + 0.5, lty = 1, col = "lightgrey")
    }
    dev.off()
  }
}

#' @title ascat.m.plotRawData
#' @description Plot tumour and germline BAF and LogR
#' @param ASCATobj an ASCAT object (e.g. data structure from ascat.loadData)
#' @param raw_LogR vector with the LogR values before correction
#' @param pch type of data points in plot
#' @param cex size of data points in plot
#' @param lim_logR y-axis limits on logR plot
#' 
#' @return Produces png files showing the logR and BAF values for tumour and germline samples
#' @author Peter Van Loo
#' @export

ascat.m.plotRawData = function(ASCATobj, raw_LogR, pch, cex, lim_logR) {
  
  print.noquote("Plotting tumor data")
  for (i in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    colls <- ifelse(ASCATobj$Germline_BAF[,i] < 0.85 & ASCATobj$Germline_BAF[,i] > 0.15, "red", "grey")
    # set point colours to show SNP germline genotype
    png(filename = paste(ASCATobj$samples[i],".tumour.png",sep=""), width = 2000, height = 1250, res = 200)
    par(mar = c(0.5,5,5,0.5), mfrow = c(3,1), cex = 0.4, cex.main=3, cex.axis = 2,
        pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000,".",20))
    plot(c(1,dim(ASCATobj$Tumor_LogR)[1]), c(-lim_logR ,lim_logR ),
         type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", tumor data, raw LogR", sep = ""), 
         xlab = "", ylab = "")
    points(ASCATobj$Tumor_LogR[,i],col="grey")
    #points(ASCATobj$Tumor_LogR[,i],col=rainbow(24)[ASCATobj$SNPpos$Chr])
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = ASCATobj$ch[[j]];
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,2,ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    plot(c(1,dim(ASCATobj$Tumor_LogR)[1]), c(-lim_logR ,lim_logR ),
         type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", tumor data, corrected LogR", sep = ""), 
         xlab = "", ylab = "")
    points(ASCATobj$Tumor_LogR[,i],col="grey")
    #points(ASCATobj$Tumor_LogR[,i],col=rainbow(24)[ASCATobj$SNPpos$Chr])
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = ASCATobj$ch[[j]];
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,2,ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    plot(c(1,dim(ASCATobj$Tumor_BAF)[1]), c(0,1), type = "n", xaxt = "n",
         main = paste(ASCATobj$samples[i], ", tumor data, BAF", sep = ""), xlab = "", ylab = "")
    points(ASCATobj$Tumor_BAF[,i], col=colls, pch = pch, cex = cex)
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk = ASCATobj$ch[[j]];
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    dev.off()
  }
  
  if(!is.null(ASCATobj$Germline_LogR)) {
    print.noquote("Plotting germline data")
    for (i in 1:dim(ASCATobj$Germline_LogR)[2]) {
      png(filename = paste(ASCATobj$samples[i],".germline.png",sep=""), width = 2000, height = 750, res = 200)
      par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2,
          pch = ifelse(dim(ASCATobj$Tumor_LogR)[1]>100000,".",20))
      plot(c(1,dim(ASCATobj$Germline_LogR)[1]), c(-1,1), type = "n", xaxt = "n",
           main = paste(ASCATobj$samples[i], ", germline data, LogR", sep = ""), xlab = "", ylab = "")
      points(ASCATobj$Germline_LogR[,i],col="red")
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk = ASCATobj$ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,2,ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      plot(c(1,dim(ASCATobj$Germline_BAF)[1]), c(0,1), type = "n", xaxt = "n",
           main = paste(ASCATobj$samples[i], ", germline data, BAF", sep = ""), xlab = "", ylab = "")
      points(ASCATobj$Germline_BAF[,i], col = colls, pch = pch, cex = cex)
      abline(v=0.5,lty=1,col="lightgrey")
      chrk_tot_len = 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk = ASCATobj$ch[[j]];
        chrk_tot_len_prev = chrk_tot_len
        chrk_tot_len = chrk_tot_len + length(chrk)
        vpos = chrk_tot_len;
        tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
        text(tpos,1,ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v=vpos+0.5,lty=1,col="lightgrey")
      }
      dev.off()
    }
  }
}

#' Plot BAF and logR profiles with ggplot
#'
#' @param dt data.frame with methylation info
#' @param outfile character srting with output pdf filename
#' Saves a pdf w/ methylation rate distribution, biases at polymorphic and 
#' non-polymorphic CG/CCGG and coverage distribution 
#' @author Elizabeth Larose Cadieux
plot_BAF_and_LogR <- function (dt, outfile, downsample=1E5) {
  
  # Only plot heterozygous SNPs
  dt_sample <- dt[dt$type=="Heterozygous",]
  sampled<-sample(1:nrow(dt_sample),size = 1E4)
  dt_sample <- dt_sample[sort(sampled),];rm(sampled)
  
  # order the dt
  dt_sample[dt_sample$chrom == "X", c("chrom")] <- 23
  dt_sample <- as.data.frame(data.table(dt_sample)[order(as.numeric(chrom), POS)])
  dt_sample[dt_sample$chrom == 23, c("chrom")] <- "X"
  dt_sample$BAFloci <- 1:length(dt_sample[,1])
  
  # get chromosome limits and labels
  SNPpos = dt_sample[,c("chrom","POS","BAFloci")]
  labels_pos <- SNPpos %>% group_by(chrom) %>% summarise_at(vars(BAFloci), 
                funs(0.5*(min(., na.rm=TRUE)+max(., na.rm=TRUE))))
  labels_pos <- unlist(unname(c(labels_pos[,2])))
  lines_pos <- data.frame(SNPpos %>% group_by(chrom) %>% summarise_at(vars(BAFloci), 
                          funs(max(., na.rm=TRUE))), stringsAsFactors = F)
  rm(SNPpos)
  
  p_LogR <- ggplot(dt_sample, aes(x=BAFloci,y=LogR_t_corr,color = flag)) + theme_minimal() +
    geom_point(shape = 20, alpha = 0.1, size = 2) + 
    annotate("text",x=labels_pos,y=rep(2.1, 23),label=factor(c(1:22,"X"), levels=c(1:22,"X")),hjust=0.5) + 
    scale_y_continuous("LogR",breaks = seq(-2.25, 2.25, by = 0.50), limits = c(-2.25, 2.25)) + ggtitle("LogR") + 
    scale_x_continuous("SNP loci", minor_breaks = lines_pos$BAFloci, breaks = lines_pos$BAFloci, labels = NULL) +
    scale_color_manual(name = "SNP\nflag", values=c("CCGG"="red","CG"="orange3","neither"="cornflowerblue"))+
    guides(color = guide_legend(override.aes = list(size=10)))

  p_BAF <- ggplot(dt_sample, aes(x=BAFloci,y=rBAF,color = flag)) +  theme_minimal() +
           geom_point(shape = 20, alpha = 0.1, size=2) + 
    annotate("text",x=labels_pos,y=rep(1.1, 23),label=factor(c(1:22,"X"), levels=c(1:22,"X")),hjust=0.5)  + ylim(0,1.1)+
    scale_color_manual(name = "SNP\nflag", values=c("CCGG"="red","CG"="orange3","neither"="cornflowerblue"))+
    scale_x_continuous("SNP loci", minor_breaks = lines_pos$BAFloci, breaks = lines_pos$BAFloci, labels = NULL) +
    ggtitle("BAF") + theme(legend.title = element_text(hjust = 0.5)) + guides(color = guide_legend(override.aes = list(size=10)))
  
  d_BAF_n <-ggplot(dt_sample, aes(x=BAF_n,y=..count..,color = flag, fill=flag)) + geom_density(alpha=0.25) +
    scale_color_manual(name = "SNP\nflag", values=c("CCGG"="red","CG"="orange3","neither"="cornflowerblue")) + 
    scale_fill_manual(name = "SNP\nflag", values=c("CCGG"="red","CG"="orange3","neither"="cornflowerblue")) + 
    theme_classic()
  
  h_BAF_n <-ggplot(dt_sample, aes(x=BAF_n,y=..count..,color = flag, fill=flag)) + geom_histogram(bins=100) +
    scale_color_manual(name = "SNP\nflag", values=c("CCGG"="red","CG"="orange3","neither"="cornflowerblue")) + 
    scale_fill_manual(name = "SNP\nflag", values=c("CCGG"="red","CG"="orange3","neither"="cornflowerblue")) + theme_classic()
  
  gglist1 <- list(p_LogR)
  gglist2 <- list(p_BAF)
  my_layout <- rbind(c(1), c(2))
  gd <- arrangeGrob(grobs = c(gglist1, gglist2), layout_matrix = my_layout)
  Sys.sleep(1)
  ggsave(file = outfile,plot = gd, device = "pdf", width=10.7, height = 6, units = "in")
  
}

#' Plot SNP data summary and QC
#' 
#' \code{plot_SNP_info} plots SNP QC
#'
#' @param dt data.table with SNP info
#' @param outfile character srting with output pdf filename
#' 
#' @return pdf
#' @author Elizabeth Larose Cadieux
plot_SNP_info <- function (dt, outfile, min) {
    
  # Total INFORMATIVE counts at SNPs
  # p1 <- ggplot(dt, aes(x=factor(1),y=total_counts), fill="black", color="black") +
  #   scale_y_log10(limits = c(min,max(dt$total_counts)+0.5*max(dt$total_counts)),
  #                 labels=scales::trans_format('log10',scales::math_format(10^.x))) +
  #   geom_violin(position = position_dodge(1), alpha = 0.25) + theme_minimal() +
  #   ylab("log(Informative coverage)") + ggtitle("A. Informative coverage") +
  #   theme(legend.title = element_text(hjust = 0.5),
  #   axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  # Total depth at SNPs (including non-informative counts)
  p2 <- ggplot(dt, aes(x=factor(1),y=total_depth), fill="black", color="black") +
        scale_y_log10(limits = c(min,max(dt$total_depth)+0.5*max(dt$total_depth)),
                      labels=scales::trans_format('log10',scales::math_format(10^.x))) +
        geom_violin(position = position_dodge(1), alpha = 0.25) + theme_minimal() +
        ylab("log(SNP coverage)") + ggtitle("A. Total SNP coverage") +
        theme(legend.title = element_text(hjust = 0.5), axis.title.x = element_blank(),
              axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  # Compare logR and logR corrected 
  #p3 <- ggplot(dt, aes(y=..count..))+
  #      ggtitle("C.")+ylab("Number of SNPs")+xlab("LogR")+
  #      theme_classic()+coord_cartesian(xlim=c(-2.5, 2.5)) +
  #      geom_histogram(aes(x=LogR_t, color="raw", fill="raw"),
  #                     binwidth = 0.05,alpha=0.05, position = "stack", size=0.75) +
  #      geom_histogram(aes(x=LogR_t_corr, fill="corrected",color="corrected"),
  #                     binwidth = 0.05, alpha=0.05, position = "stack",size=0.65) +
  #      scale_color_manual(name = "LogR legend",
  #                         values=c("raw"="hotpink","corrected"="chartreuse3"))+
  #      theme(legend.title = element_text(hjust = 0.5, size = 10, face = "bold"),
  #            legend.text = element_text(size = 10))+
  #      scale_fill_manual(name = "LogR legend",
  #                        values=c("raw"="hotpink","corrected"="chartreuse3"))
  
  # Look at nuumber of CCGG,CG and other SNPs
  tmp <- data.frame(table(dt$flag))
  colnames(tmp) <- c("flag", "frequency")
  tmp$flag <- factor(tmp$flag,levels = c("CCGG","CG","neither"), ordered = T)
  tmp$labels <- paste0(round(tmp$frequency/sum(tmp$frequency)*100, digits=0),"%")

  p4 <- ggplot(tmp, aes(x=factor(1), y=frequency, fill=flag, color=flag,
               label=ifelse(round(frequency/sum(frequency)*100, digits=0)<1,NA,labels)))+
    theme_classic() + ggtitle("C. SNPs context") +
    geom_bar(width = 1, alpha=0.25, stat="identity") +
    geom_text(size = 3, position = position_stack(vjust = 0.5), color="black") +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.line = element_blank(),legend.text=element_text(size=8),
          axis.ticks = element_blank()) +
    scale_color_manual(name = "SNP flag",
                       values=c("CCGG"="darkblue","CG"="darkred","neither"="slategray"))+
    scale_fill_manual(name = "SNP flag",
                      values=c("CCGG"="darkblue","CG"="darkred","neither"="slategray"))+
    coord_polar("y") 
  
  # Look at nuumber of CCGG,CG and other SNPs
  tmp <- data.frame(table(dt$type))
  colnames(tmp) <- c("type", "frequency")
  tmp$type <- factor(tmp$type,levels = c("Homozygous" ,"Heterozygous"), ordered = TRUE)
  tmp$labels <- paste0(round(tmp$frequency/sum(tmp$frequency)*100, digits=0),"%")
  
  p5 <- ggplot(tmp, aes(x=factor(1), y=frequency, fill=type, color=type,
               label=ifelse(round(frequency/sum(frequency)*100, digits=0)<1,NA,labels)))+
        theme_classic() + ggtitle(paste0("E. Germline genotype")) +
        geom_bar(width = 1, alpha=0.25, stat="identity") +
        geom_text(size = 3, position = position_stack(vjust = 0.5), color="black") +
        theme(axis.text = element_blank(), axis.title = element_blank(),
              axis.line = element_blank(),legend.text=element_text(size=8),
              axis.ticks = element_blank()) +
        scale_color_manual(name = "Germline\nSNP genotype",
              values = c("Homozygous" = "orange3", "Heterozygous" = "mediumpurple")) +
        scale_fill_manual(name = "Germline\nSNP genotype",
              values = c("Homozygous" = "orange3", "Heterozygous" = "mediumpurple")) +
        coord_polar("y")
    
  # Look at number of homo/het SNPs at CCGG,CG and other SNPs
  tmp <- data.frame(table(dt[, .SD, .SDcols=c("flag", "type")]))
  colnames(tmp) <- c("SNP Context", "Germline Genotype", "Number of SNPs")
  tmp$"SNP Context" <- factor(tmp$"SNP Context",
                              levels = c("CCGG","CG","neither","combined"), ordered = T)
  tmp$"Percentage" <- paste0(round(tmp$"Number of SNPs"/sum(tmp$"Number of SNPs")*100,
                             digits=2))
  
  tmp2 <- data.frame(table(dt[, .SD, .SDcols=c("flag")]))
  tmp2$"Germline Genotype" <- NA ; tmp2$"Germline Genotype" <- "combined" 
  colnames(tmp2) <- c("SNP Context", "Number of SNPs", "Germline Genotype" )
  tmp2 <- tmp2[,c(1,3,2)]
  tmp2$"SNP Context" <- factor(tmp2$"SNP Context",
                               levels = c("CCGG","CG","neither","combined"), ordered = T)
  tmp2$"Percentage" <- paste0(round(tmp2$"Number of SNPs"/sum(tmp2$"Number of SNPs")*100,
                              digits=2))
  
  tmp3 <- data.frame(table(dt[, .SD, .SDcols=c("type")]))
  tmp3$"SNP Context" <- NA ; tmp3$"SNP Context" <- c("combined") 
  colnames(tmp3) <- c("Germline Genotype","Number of SNPs", "SNP Context")
  tmp3 <- tmp3[,c(3,1,2)]
  tmp3$"SNP Context" <- factor(tmp3$"SNP Context",
                               levels = c("CCGG","CG","neither","combined"), ordered = T)
  tmp3$"Percentage" <- paste0(round(tmp3$"Number of SNPs"/sum(tmp3$"Number of SNPs")*100,
                              digits=2))
  
  tmp = rbind.data.frame(tmp, tmp2, tmp3, stringsAsFactors = F)
  
  p9 <- tableGrob(tmp,rows=NULL)
  
  gglist1 <- list(p2)
  gglist2 <- list(p4,p5,p9)
  my_layout <- rbind(c(1,2,4,4),c(1,3,4,4))
  gd <- arrangeGrob(grobs = c(gglist1, gglist2), layout_matrix = my_layout)
  Sys.sleep(1)
  ggsave(file = outfile,plot = gd, device = "pdf", width=11.5, height = 4, units = "in")
  
  # #clean-up
  # rm(p,p2,p3,p4,gglist1,gglist2,gd,min,my_layout,p5,dt)
  
}

#' Plot plots SNP QC
#'
#' @param dt data.table with SNP info
#' @param outfile character srting with output pdf filename
#' 
#' @return pdf
#' @author Elizabeth Larose Cadieux
plot_normal_SNP_info <- function (dt, outfile, min) {
  tmp <- dt[BAF>=0.15 & BAF <= 0.85 & !is.na(BAF),]
  tmp2 <- table(cut(tmp$BAF, breaks = (0.85-0.15)/0.01))
  tmp2 <- unname(tmp2[which.max(tmp2)])*9/10
  p4 <- ggplot(data = tmp, aes(x=BAF,y=..count..,color = type,fill = type)) +
    geom_histogram(binwidth = 0.01, alpha = 0.25) + theme_minimal() +
    scale_color_manual(name = "", values = c("Homozygous" = "orange3", "Heterozygous" = "mediumpurple")) +
    scale_fill_manual(name = "", values = c("Homozygous" = "orange3", "Heterozygous" = "mediumpurple")) +
    geom_vline(data = tmp, aes(xintercept=mean(BAF)), color="red", linetype="dashed") +
    annotate("text",x=0.51,y=tmp2,
             label=paste("mean= ",as.character(round(mean(tmp$BAF, na.rm=TRUE), digits = 3)),
                         "\nB-allele bias= ",100*round(mean(0.5 - tmp$BAF, na.rm=TRUE), digits = 3), 
                         "%", sep=""),
             hjust=0, size=3) + ylim(c(0,nrow(tmp[BAF < 0.505 & BAF > 0.49,1]))) +
    ggtitle("A.") + theme(legend.position="none")
  
  p5 <- ggplot(dt, aes(x=type,y=total_counts,color = type,fill = type)) +
    scale_y_log10(limits = c(min,max(dt$total_counts)+0.5*max(dt$total_counts)),
                  labels=scales::trans_format('log10',scales::math_format(10^.x))) +
    geom_violin(position = position_dodge(1), alpha = 0.25) +
    theme_minimal() +
    ylab("log(Informative SNP coverage)") + xlab("Germline SNP genotype") +
    scale_color_manual(name = "Germline\nSNP genotype",
                       values = c("Homozygous" = "orange3", "Heterozygous" = "mediumpurple")) +
    scale_fill_manual(name = "Germline\nSNP genotype",
                      values = c("Homozygous" = "orange3", "Heterozygous" = "mediumpurple")) +
    ggtitle("B.") + theme(legend.title = element_text(hjust = 0.5)) +
    annotate("text",x=c("Homozygous", "Heterozygous"),
             y=rep(dt[, max(total_counts)+0.2*max(total_counts)], 2),
             label=paste(c(dt[, sum(type == "Homozygous", na.rm=TRUE)],
                           dt[, sum(type == "Heterozygous", na.rm=TRUE)]), "\nSNPs",sep=" "), vjust=0.5)
  
  dt[, SNP := paste(ref, alt, sep="")]
  dt[, SNP := factor(ifelse(SNP=="TG", "AC", ifelse(SNP=="TC", "AG",
                     ifelse(SNP=="TA", "AT", ifelse(SNP=="GA", "CT",
                     ifelse(SNP=="GC", "CG", ifelse(SNP=="GT", "CA", SNP)))))))]
  p6 <- ggplot(dt, aes(x=SNP,y=total_counts,color = SNP,fill = SNP)) +
        scale_y_log10(limits = c(min,max(dt$total_counts)+0.5*max(dt$total_counts)),
                      labels=scales::trans_format('log10',scales::math_format(10^.x))) +
    geom_violin(position = position_dodge(1), alpha = 0.25) + theme_minimal() +
    ylab("log(Informative SNP coverage)") + xlab("Germline SNP genotype") +
    scale_color_manual(name = "Germline\nSNP genotype",
                       values = c("AC"="orange3","AG"="mediumpurple","AT"="firebrick",
                                  "CA"="lightblue","CG"="darkblue","CT"="grey")) +
    scale_fill_manual(name = "Germline\nSNP genotype",
                      values = c("AC"="orange3","AG"="mediumpurple","AT"="firebrick",
                                 "CA"="lightblue","CG"="darkblue","CT"="grey")) +
    theme(legend.title = element_text(hjust = 0.5)) 
  
  tmp <- data.frame(table(dt$SNP))
  colnames(tmp)[1:2] <- c("SNP", "SNP_freq")
  tmp$anno_start <- c(0, cumsum(tmp$SNP_freq[1:(nrow(tmp)-1)]))
  tmp$anno_end <- cumsum(tmp$SNP_freq)
  tmp$mid_point <- (tmp$anno_end-tmp$anno_start) / 2 + tmp$anno_start
  tmp$labels <- paste0(round(tmp$SNP_freq/sum(tmp$SNP_freq)*100, digits=0),"%")
  
  p7 <- ggplot(tmp, aes(x=factor(1), y=SNP_freq, fill=SNP, color=SNP,
               label=ifelse(round(SNP_freq/sum(SNP_freq)*100, digits=0)<1,NA,labels)))+
    theme_classic() + geom_bar(width = 1, alpha=0.25,stat = "identity") +  
    ggtitle(paste0("C. ", sum(tmp$SNP_freq)," SNPs")) +
    geom_text(size = 3, position = position_stack(vjust = 0.5), color="black") +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.line = element_blank(),legend.text=element_text(size=8),
          axis.ticks = element_blank()) +
    scale_color_manual(name = "Germline\nSNP genotype",
                       values = c("AC"="orange3","AG"="mediumpurple","AT"="firebrick",
                                  "CA"="lightblue","CG"="darkblue","CT"="grey")) +
    scale_fill_manual(name = "Germline\nSNP genotype",
                      values = c("AC"="orange3","AG"="mediumpurple","AT"="firebrick",
                                 "CA"="lightblue","CG"="darkblue","CT"="grey")) +
    coord_polar("y") + theme(legend.position="none")
  
  gglist1 <- list(p4,p5)
  gglist2 <- list(p7,p6)
  my_layout <- rbind(c(1,1,1,2,2,2,2),c(3,3,rep(4,5)))
  gd <- arrangeGrob(grobs = c(gglist1,gglist2), layout_matrix = my_layout)
  Sys.sleep(1)
  ggsave(file = outfile,plot = gd, device = "pdf", width=7, height = 6.2, units = "in")
  
  # #clean-up
  # rm(p4,p5,p6,gglist2,gglist1,gd,my_layout,tmp)
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
