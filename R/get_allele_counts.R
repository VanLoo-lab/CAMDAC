##  CAMDAC get allele counts at SNPs and CpGs from bisulfite sequencing BAM file

#' Compile allele counts at SNPs and at CpGs for bisulfite sequencing data
#' \code{get_allele_counts}
#'
#' @param i Integer loop index. The function must be run with all values from 1 to 25, each containing
#' 1/25th of the RRBS covered genome.
#' @param patient_id Character variable containting the patient id
#' @param sample_id Character variable with the sample id
#' @param sex Character variable with the patient sex expressed as "XX" for female or "XY" for male.
#' @param bam_file Character variable with the full bam file name and path
#' @param mq Character variable or numeric containting the mapping quality treshold to be used.
#' For RRBS, set mq=0. Read mapping validity is based on read start site and nucleotides rather than mq.
#' @param path Character path variable pointing to the desired working directory.
#' This is where the output will be stored and should be constant for all CAMDAC functions.
#' Do not alter the output directory structure while running CAMDAC.
#' The function output of this function will be a sub-directory of the path variable under
#' "./Allelecounts/sample_id/". Do not change the directory structure as subsequent functions will
#' look for files in this directory.
#' @param path_to_CAMDAC Character variable containting the CAMDAC installation path  (e.g. "/path/to/CAMDAC/").
#' @param build Character variable corresponding to the reference genome used for alignment.
#' CAMDAC is compatible with "hg19", "hg38", "GRCH37","GRCH38".
#' @param n_cores Numerical value correspdonding to the number of cores for parallel processing
#' @param test Logical value indicating whether this is a quick test run with data subsampling
#' 
#' @return One .fst file including methylation info at CpGs and BAF and depth of coverage at
#' SNPs for the ith subset of RRBS loci

get_allele_counts <- function (i , patient_id, sample_id, sex, bam_file, mq=0,
                               path, path_to_CAMDAC, build=NULL, n_cores, test=FALSE){
  
  if(getOption("scipen")==0){options(scipen = 999)} 
  # important to turn scientific notation off when saving genomic coordinates to .txt files
  
  # ensure mq is parsed as numerical value
  mq <- as.numeric(mq)
  cat("Mapping treshold MQ ≥ ",mq," applied","\nBase quality treshold BQ ≥ 20 applied\n", sep = "")
  
  # Set working directory path and create results folders
  # Do not change this - subsequent functions will look for files in this directory
  path_output <- file.path(path, patient_id, "Allelecounts", sample_id)
  suppressWarnings(dir.create(path_output, recursive = TRUE))

  # Load doParrellel if running job in parrallel 
  if(n_cores > 1){
    x <- c("doParallel", "parallel")
    invisible(lapply(x, function(y) suppressWarnings(suppressMessages(library(y,
    character.only = TRUE,quietly=TRUE, warn.conflicts = FALSE)))));rm(x)
  }
  
  # Ensure bam file and bai file exists
  if(!file.exists(bam_file)){
    errorCondition(message = ".bam file not found at path provided - please provide bam file name including the full path")
  }
    
  # Detect build if not provided
  if(is.null(build)){
    tmp <- scanBamHeader(bam_file)[[1]][[1]]
    tmp <- tmp[grepl("X$",names(tmp))]
    build <- ifelse(unname(tmp)=="155270560", "hg19", "hg38")
    UCSC <- ifelse(substr(names(tmp), 1, 1)=="c", TRUE, FALSE)
  }
  cat(paste("ScanBam pileup with build ", build, sep = " "), "\n", sep = "")
  
  # Set build variables
  if(build%in%c("hg19","hg38")){UCSC=TRUE}
  if(build%in%c("GRCH37","GRCH38")){UCSC=FALSE}
  if(build=="GRCH37"){build="hg19"} # set build to to assembly version disregarging USCS/Ensembl
  if(build=="GRCH38"){build="hg38"}

  # set path to reference files
  loci_file=file.path(path_to_CAMDAC,paste0("pipeline_files/loci_files/"))
  segments_file_path=file.path(path_to_CAMDAC,"pipeline_files/segments_files/")
  
  # load segments file
  f_name = paste(segments_file_path, "segments.",build, ".",i, ".RData", sep = "")
  load(f_name)
  rm(f_name)
  
  # Ensure that spurious alignments to Y in females are removed
  if(sex=="XX"&i==25){segments_subset<-segments_subset[!as.character(seqnames(segments_subset))%in%c("chrY","Y")]}
  
  # Modify seqlevels if necessary
  if(UCSC == FALSE){ seqlevelsStyle(segments_subset) <- "Ensembl" }
  
  # Force relevant seqlevels
  if(UCSC == FALSE){chr_names = c(1:22, "X", "Y")}
  if(UCSC == TRUE){chr_names = paste0("chr",c(1:22, "X", "Y"))}
  seqlevels(segments_subset) <- chr_names
  
  # function to get SNP and CpG allele counts
  get_reads <- function(i, j = NULL, n_cores = NULL, bam_file, segments_subset, mq, test){

    # Store segements file in tmp object
    segments_subset_2 <- segments_subset
    
    # If this is a test run - downsample
    if(test==TRUE){
      set.seed(0)
      segments_subset_2 <- sort(segments_subset_2[sample(1:length(segments_subset_2), 5E3)])
      }
    
    # If running in parallel, further subset the genomic segments file (input for scanBam)
    if(!is.null(j)){
      segments_subset_2$coreID <- as.numeric(ceiling(1:length(start(segments_subset_2))/
                                            (length(start(segments_subset_2))/n_cores+1)))
      segments_subset_2 <- segments_subset_2[segments_subset_2$coreID==j,]
    }
    
    # remove metadata
    mcols(segments_subset_2) <- NULL
    
    # Set what argument to speed up scanBam
    whats <- c("qname", "rname", "strand", "pos", "qwidth", "mapq", "seq", "qual") # scanBamWhat()
    idxFile <- paste(bam_file, ".bai", sep = "")
    
    if(!file.exists(idxFile)){
      errorCondition(message = ".bam.bai file not found at path provided - 
                     please ensure the .bam.bai file is in the same directory as your 
                     bam file and has the same file name with added .bai extension")
    }
    
    # Set info that we want to obtain for each position
    # names(scanBamHeader(bam_file)[[1]][[1]])[1] == "chr1"
    bamFile <- BamFile(bam_file)
    params <- Rsamtools::ScanBamParam(which=segments_subset_2, what=whats)
    # Rsamtools::scanBamFlag(isPaired = TRUE, isFirstMateRead = TRUE)
    # which <- GRanges(seqnames = c("chrY"), ranges = IRanges(1,16569))
    
    # Run ScanBam 
    #bam <- Rsamtools::scanBam(file=bam_file, index=idxFile, param=params, flag = scanBamFlag())
    bam <- Rsamtools::scanBam(file=bamFile, param=params)
    rm(params, segments_subset_2)
    
    # Remove entries with zero read
    bam_filtered <- bam[which(sapply(bam, function(x) length(x[["qname"]])) != 0)]
    rm(bam)
      
    # Unlist bam
    df <- data.table(dplyr::bind_rows(lapply(bam_filtered, as.data.frame.list, stringsAsFactors=FALSE)))
    # duplicates created when unlisting and due to overlapping reads, to be removed in later lines
    rm(bam_filtered)
    
    # Return empty BAM file if no data
    if(nrow(df) == 0){
      df_names = c("CHR","chrom","start","end","width",
                   "POS","ref","alt","alt_counts","ref_counts",
                   "total_counts","BAF","total_depth",
                   # "other_counts","all_counts",
                   "M","UM","total_counts_m","m",
                   # allele counts breakdown may be obtained by uncommenting this line + line 441.
                   #"Af","Ar","Cf","Cr","Tf","Tr","Gf","Gr","CAr","TGf","CGf","CGr", 
                   "CCGG")
      df <- data.frame(matrix(NA, nrow=1, ncol=length(df_names)))
      names(df) <- df_names
      df <- df[-1,] # Return empty dataframe
      return(df)
    } 
        
    # Format pileup df for transformation into GRanges Object
    df$end <- df$pos + (df$qwidth - 1)
    colnames(df)[c(which(colnames(df) == "rname"),which(colnames(df) == "pos"))] <- c("chrom", "start")
    rownames(df) <- c(1:nrow(df))
  
    cols <- c("qname", "chrom", "strand", "start", "end", "mapq", "seq", "qual")
    Sys.sleep(1)
        
    df_bam <- df[,cols,with=FALSE] 
    rm(df)
        
    # The pileup unlisting produces duplicate rows where more than one loci overlaps with the same molecule
    x <- do.call(paste, c(df_bam[,cols,with=FALSE], sep="-")) 
    df_bam_dedup <- df_bam[!duplicated(x),]
    rm(df_bam, cols,x)
 
    # Create Granges - Each row represents a single read 
    gr_bam <- makeGRangesFromDataFrame(df_bam_dedup, keep.extra.columns = TRUE, na.rm=TRUE)
    rm(df_bam_dedup)
        
    # Load reference CpG/SNP loci file
    load(file.path(loci_file, paste0("loci.", build, ".",i, ".RData")))

    # Change seqnames if necessary
    if(UCSC == FALSE){seqlevelsStyle(loci_subset) <- "Ensembl"}
    
    # Extract info per CpG/SNP loci
    overlaps <- IRanges::mergeByOverlaps(gr_bam, loci_subset)
    rm(gr_bam, loci_subset)
        
    df_pileup <- data.table(qname=as.character(overlaps$qname), strand=as.character(strand(overlaps[,"gr_bam"])),
                            chrom=as.character(seqnames(overlaps[,"gr_bam"])),read.start=start(overlaps[,"gr_bam"]), 
                            read.end=end(overlaps[,"gr_bam"]),POS=overlaps$POS, width = width(overlaps[,"loci_subset"]),
                            start=start(overlaps[,"loci_subset"]),end=end(overlaps[,"loci_subset"]),ref=overlaps$ref,
                            alt=overlaps$alt,seq=overlaps$seq,qual=overlaps$qual,mq=overlaps$mapq,
                            stringsAsFactors = FALSE)
    rm(overlaps)
        
    # Make sure overlapping mates at a given position are only counted once for paired-end data
    # The qname and CpG/SNP position will be the same for overlapping RRBS paired-end reads
    cols<-c("qname","chrom","POS","width","start","ref","alt")
    x <- do.call(paste, c(df_pileup[,cols,with=FALSE], sep="-")) 
    df_pileup_dedup <- df_pileup[!duplicated(x),]
    rm(df_pileup, cols, x)
    
    # Flag 5'CCGG and rm 3' CCGG but keep any SNPs. 
    df_pileup_dedup[, CCGG := (ifelse(width != 4 , "non-CCGG",
                               ifelse(((strand == "+" & read.start < end & read.start >= start) | 
                                       (strand == "-" & read.end >= start & read.end <= end)), "5pCCGG",
                               ifelse(((strand == "-" & read.start <= end & read.start >= start) |
                                       (strand == "+" & read.end >= start & read.end <= end)), "3pCCGG",
                               ifelse(!is.na(POS), "CCGG-SNP", "CCGG-SNV")))))]
  
    # Assuming complete digestion, non-3p/5p CCGG can occur due to: 
    # - a germline/tumour mutation destroying CCGG motif (CCGG-SNV) or
    # - a SNP allele destroying the CCGG motif (CCGG-SNP) 
    # Problems with NuGEN adapter trimming could also lead to within read CCGGs.
     
    # Remove mid-read CCGGs if they don't overlap with a 1000 genome SNPs. 
    df_pileup_dedup[, flag := ifelse(CCGG=="CCGG-SNV",FALSE,TRUE)]
    df_pileup_dedup <- df_pileup_dedup[(flag),] 
    df_pileup_dedup[, flag := NULL] 
    
    # Flag SNPs at CCGG-destroying SNPs 
    df_pileup_dedup[, flag := (width==4 & !is.na(POS))] 
    df_pileup_dedup[, CCGG_SNP := (ifelse(flag == FALSE, "n/a",
                                   ifelse(CCGG=="5pCCGG" & # Flag 5' SNPs inside of read 
                                         ((strand=="+" & POS!=start) | (strand=="-" & POS!=end)), "yes",
                                   ifelse(CCGG=="3pCCGG" & # Flag 3' SNPs inside of read 
                                         ((strand=="+" & POS==start) | (strand=="-" & POS==end)), "yes", 
                                   ifelse(CCGG=="CCGG-SNP","yes","no")))))]
    df_pileup_dedup[, flag := NULL]
    
    # Extract 5'CCGG alleles
    df_pileup_dedup[, flag := (CCGG=="5pCCGG")]
    df_pileup_dedup[, tmp := (start-read.start)]
    df_pileup_dedup[, alleles.CCGG := ifelse(flag == FALSE, NA, 
                                       paste(substr(as.character(seq),tmp+1,tmp+4),strand,sep=""))]
    df_pileup_dedup[, flag := NULL]
    
    # Annotate allowed 5' alleles CGG+, TGG+, CCG- or CCA- 
    z <- grep("(^CC([AG]{1})([ACGNT]{0,1})([-])$)|(([ACGNT]{0,1})([CT]{1})GG([+]))",df_pileup_dedup$alleles.CCGG)
        
    # Remove reads where 5' CCGG isn'TRUE match [CT]CGG+, [CT]TGG+, CCG[GA]- or CCA[GA]- 
    # These are usually due to poor alignment (i.e. tend to be C rich )
    df_pileup_dedup[, flag := FALSE]
    df_pileup_dedup[z, flag := TRUE] ; rm(z)
    
    # Find all misaligned reads
    qnames <- unique(df_pileup_dedup[,(qname[flag==FALSE & CCGG=="5pCCGG"])])
    
    # Keep reads which have an alternative overlapping 5'CCGG w/ suitable alleles 
    tmp <- unique(df_pileup_dedup[,(qname[qname%in%qnames & flag==TRUE])])
    qnames <- qnames[!qnames%in%tmp];rm(tmp)
    
    # Remove misaligned reads
    df_pileup_dedup<-df_pileup_dedup[!qname%in%qnames,];rm(qnames)
    
    # Rename remaining non-5p CCGGs which are due to SNP / SNV
    flag2 <- df_pileup_dedup[, (flag==FALSE & CCGG=="5pCCGG")]
    df_pileup_dedup[flag2 & is.na(POS), CCGG := "CCGG-SNV"] 
    df_pileup_dedup[flag2 & !is.na(POS), CCGG := "CCGG-SNP"] 
    # Note. SNPs/SNVs are unlikely to form a valid 5p CCGG.
    # To do so, it would have to turn a motif overlappingto the ref CCGG to another CCGG 
    # (i.e. -CCGGCGG- to -CCGCCGG-)
    
    # clean-up
    df_pileup_dedup[, flag := NULL] ; rm(flag2)
    
    # # print diganostics
    # cat("Misaligned reads removed", "Adapter contamination removed!", sep = "\n", file = file_stdout_2, append = TRUE)
    
    # Extract CpG and SNP alleles
    df_pileup_dedup[, offset := ifelse(width == 4, 2, 1)]
    df_pileup_dedup[, alleles.dinucs :=  ifelse(width>=2&CCGG!="3pCCGG", 
                             paste(substr(seq, tmp+offset, tmp+offset+1),strand, sep=""), NA)]
    df_pileup_dedup[, tmp2 := (POS-read.start+1)]
    df_pileup_dedup[, alleles.SNP := ifelse(!is.na(POS)&CCGG_SNP!="no", paste(substr(seq,tmp2,tmp2),strand, sep=""), NA)]
    
    # Extract CpG and SNP quality from the base quality string
    df_pileup_dedup[, qual.dinucs := ifelse(width>=2&CCGG!="3pCCGG",substr(qual,tmp+offset, tmp+offset+1), NA)]
    df_pileup_dedup[, qual.SNP := ifelse(!is.na(POS)&CCGG_SNP!="no",substr(qual,tmp2,tmp2),NA)]
      
    # Concatenate results
    cols <- c("chrom","strand", "start", "end", "width","CCGG", "alleles.CCGG", "POS", "ref", "alt",
              "alleles.dinucs","alleles.SNP","qual.dinucs","qual.SNP","mq")
    assign("df_pileup_alleles",df_pileup_dedup[,cols,with=FALSE]);rm(df_pileup_dedup,cols)
    
    # Flag SNPs and CpGs where BQ < 20 
    df_pileup_alleles[, flag1 := width >= 2 & grepl(pattern = "([5-9A-K:-@])([5-9A-K:-@])", qual.dinucs)]
    df_pileup_alleles[, flag2 := !is.na(POS) & grepl(pattern = "([5-9A-K:-@])", qual.SNP)]
    df_pileup_alleles[, flag := (flag1==TRUE | flag2==TRUE)] # CpG or SNP has BQ > 20
   
    # Remove positions where neither CG and/or SNP passes the BQ tresholds
    if(sum(df_pileup_alleles$flag)!=nrow(df_pileup_alleles)){
      df_pileup_alleles <- df_pileup_alleles[df_pileup_alleles$flag == TRUE,] }
    df_pileup_alleles[, flag := NULL]
        
    # Remove dinucs at CpGs where only SNPs passes qual treshold, but keep SNP loci for BAF and logR
    df_pileup_alleles[, alleles.dinucs := ifelse(flag1==TRUE, alleles.dinucs, NA)]
    df_pileup_alleles[, qual.dinucs := ifelse(flag2==TRUE, qual.dinucs, NA)]
    df_pileup_alleles[, c("flag1", "flag2") := NULL]
    
    # remove reads where MQ < mq treshold (recommeded is 0)
    if(mq > 0){df_pileup_alleles <- df_pileup_alleles[df_pileup_alleles$mq >= mq,]}
  
    # # print diganostics
    # cat("Mapping treshold MQ ≥ ",mq," applied","\nBase quality treshold BQ ≥ 20 applied\n", sep = "", 
    #     file = file_stdout_2, append = TRUE)
    
    # Assign ref or alt count for each read
    df_pileup_alleles[, Af := ifelse(alleles.SNP == "A+",1,0)]
    df_pileup_alleles[, Ar := ifelse(alleles.SNP == "A-",1,0)]
    df_pileup_alleles[, Cf := ifelse(alleles.SNP == "C+",1,0)]
    df_pileup_alleles[, Cr := ifelse(alleles.SNP == "C-",1,0)]
    df_pileup_alleles[, Gf := ifelse(alleles.SNP == "G+",1,0)]
    df_pileup_alleles[, Gr := ifelse(alleles.SNP == "G-",1,0)]
    df_pileup_alleles[, Tf := ifelse(alleles.SNP == "T+",1,0)]
    df_pileup_alleles[, Tr := ifelse(alleles.SNP == "T-",1,0)]
  
    # Assign ref or alt count for each read
    df_pileup_alleles[, CGr := ifelse(alleles.dinucs == "CG-",1,0)]
    df_pileup_alleles[, CGf := ifelse(alleles.dinucs == "CG+",1,0)]
    df_pileup_alleles[, CAr := ifelse(alleles.dinucs == "CA-",1,0)]
    df_pileup_alleles[, TGf := ifelse(alleles.dinucs == "TG+",1,0)]
  
    # Assign fragments breakpoint
    df_pileup_alleles[, CCGG := ifelse(CCGG == "5pCCGG",1,0)]
    
    # Calculate total ref or alt count for each position
    df_summary <- df_pileup_alleles[,.("Af" = sum(Af, na.rm = TRUE), "Ar" = sum(Ar, na.rm = TRUE),
                                       "Cf" = sum(Cf, na.rm = TRUE), "Cr" = sum(Cr, na.rm = TRUE),
                                       "Gf" = sum(Gf, na.rm = TRUE), "Gr" = sum(Gr, na.rm = TRUE),
                                       "Tf" = sum(Tf, na.rm = TRUE), "Tr" = sum(Tr, na.rm = TRUE),
                                       "CGf" = sum(CGf, na.rm = TRUE), "CGr" = sum(CGr, na.rm = TRUE),
                                       "TGf" = sum(TGf, na.rm = TRUE), "CAr" = sum(CAr, na.rm = TRUE),
                                       "total_depth" = length(Af), "CCGG" = sum(CCGG, na.rm = TRUE)),
                                       #"mq"=paste0(mq,collapse=",")),
                                   keyby = .(chrom, start, end, width, POS, ref, alt)]
    rm(df_pileup_alleles)
    
    # Remove positions with only low mapping quality reads
    # df_summary$flag <- df_summary$mean_mq>=13
    # if the mean mapping quality < 13 remove the entire Msp1 frgment 
    # (i.e. prob that read is mapped incorrectly > 5% )
    
    # # print diganostics
    # cat("Allele pileup: complete!", sep = "\n", file = file_stdout_2, append = TRUE)
        
    #Get SNP in ref:alt format 
    df_summary[, SNP := ifelse(!is.na(ref), paste(ref, alt, sep = ""), NA)]
    df_summary[, offset := ifelse(width == 4, 1, 0)]
  
    #Calculate the ref counts
    df_summary[, ref_counts := ifelse(is.na(ref), NA,
                               ifelse(grepl(pattern="AC",SNP), Af + Ar,
                               ifelse(grepl(pattern="CA",SNP), Tf + Cr + Cf,
                               ifelse(grepl(pattern="AG",SNP), Af,
                               ifelse(grepl(pattern="GA",SNP), Gf,      
                               ifelse(grepl(pattern="AT",SNP), Af + Ar,
                               ifelse(grepl(pattern="TA",SNP), Tf + Tr,
                               ifelse(grepl(pattern="GT",SNP), Gf + Ar + Gr,
                               ifelse(grepl(pattern="TG",SNP), Tf + Tr,
                               ifelse(grepl(pattern="CG",SNP), Tf + Cr + Cf,
                               ifelse(grepl(pattern="GC",SNP), Gf + Ar + Gr,
                               ifelse(grepl(pattern="CT",SNP), Cr,
                               ifelse(grepl(pattern="TC",SNP), Tr,
                                    NA)))))))))))))]
    
    #We consider all bisulfite converted alleles and that C methylation is possible both in and outside CpG context
    
    #Calculate the alt counts
    df_summary[, alt_counts := ifelse(is.na(ref), NA,
                               ifelse(grepl(pattern="AC",SNP), Tf + Cr + Cf,
                               ifelse(grepl(pattern="CA",SNP), Af + Ar,
                               ifelse(grepl(pattern="AG",SNP), Gf,
                               ifelse(grepl(pattern="GA",SNP), Af,      
                               ifelse(grepl(pattern="AT",SNP), Tf + Tr,
                               ifelse(grepl(pattern="TA",SNP), Af + Ar,
                               ifelse(grepl(pattern="GT",SNP), Tf + Tr,
                               ifelse(grepl(pattern="TG",SNP), Gf + Ar + Gr,
                               ifelse(grepl(pattern="CG",SNP), Gf + Ar + Gr,
                               ifelse(grepl(pattern="GC",SNP), Tf + Cr + Cf,
                               ifelse(grepl(pattern="CT",SNP), Tr,
                               ifelse(grepl(pattern="TC",SNP), Cr, 
                                    NA)))))))))))))]          
    
    df_summary[, total_counts := ifelse(is.na(POS), NA, ref_counts+alt_counts)]
    
    # Counts methylated for allele-specific methylation i.e. excluding CG destroying SNPs.
    df_summary[, M := ifelse(width == 1, NA,
                      ifelse(is.na(SNP), CGf + CGr,
                      ifelse(start+offset == POS & grepl(pattern = "([CT][TC])",SNP), CGr,
                      ifelse(end-offset == POS &grepl(pattern = "([GA][AG])",SNP), CGf,
                             CGf + CGr))))]
    
    # Get the total methylation counts --> CG-destroying SNP allele does not contribute. 
    # Will be NA if the CG-destroying allele at CG/CCGG is the only one present
    df_summary[, total_counts_m := ifelse(width == 1, NA,
                                   ifelse(is.na(SNP), TGf + CAr + M,
                                   ifelse(start+offset == POS & grepl(pattern = "([CT][TC])", SNP), CAr + M,
                                   ifelse(end-offset == POS & grepl(pattern = "([GA][AG])", SNP), TGf + M, 
                                          TGf + CAr + M))))]
    
    # Get the counts unmethylated
    df_summary[, UM := total_counts_m-M]
    df_summary[, offset := NULL]
    
    # Where there is insufficient evidence that a CG motif is present or what its methylation status is, set M and UM to NA
    df_summary[total_counts_m <= 3 & width >= 2, c("M","UM","total_counts_m") := NA ]
    
    # Calculate the all counts including ones indistinguishable between alt and ref and /or meth and unmethyl
    df_summary[, all_counts := ifelse(is.na(ref), TGf+CAr+CGf+CGr,
                               ifelse(grepl(pattern="([GA][AG])",SNP),Af+Gf+Ar+Gr,
                               ifelse(grepl(pattern="([CT][TC])",SNP),Cr+Tr+Tf+Cf,
                                      total_counts)))]
    
    # Over-conversion of methylated CGs has not been accounted for and should be in later versions.
    # Under-conversion is to be ignored. At SNPs -> less than 1 in 499 reads affected.
  
    # total_counts == expected nucletoides distinguishable ref and alt nucs only
    # total_counts_m == expected dinucletoides distinguishable between meth / unmeth (excluding CG destroying SNPs)
    # all_counts == expected dinucletoides/nucleotide including indistinguishable ref and alt and or meth/unmeth bases 
    # total_depth == all reads
    # other_counts == unexpected nucleotides
    
    #Remove rows with no allowed counts
    flag <- df_summary[, (is.na(total_counts_m)&is.na(total_counts))]
    if(sum(flag)>0){df_summary <- df_summary[flag==FALSE,]};rm(flag)
    
    # Remove positions with unexpected bases at SNPs (i.e. SNV)
    df_summary[, other_counts :=  total_depth - all_counts]
    flag <- df_summary[, (ifelse(other_counts <= 0.05*total_depth, TRUE, FALSE))]
      if(!is.na(table(flag)["FALSE"])) {
          warning(paste0("Unexpected nucleotides represent > 5% of the alleles at ", sum(table(flag)["FALSE"]), 
                         " positions ( ", round(sum(table(flag)["FALSE"])*100/length(flag), digits = 2),
                         " %) - positions removed.", sep = ""))
          df_summary<- df_summary[flag,]
      };rm(flag)

    # Select and order cols    
    colls <- c("chrom","start","end","width","POS",
               "ref","alt","alt_counts","ref_counts","total_counts","total_depth",
               #"other_counts","all_counts",
               "M","UM","total_counts_m",
              # allele counts breakdown may be obtained by uncommenting this line + line 515.
               #"Af","Ar","Cf","Cr","Tf","Tr","Gf","Gr","CAr","TGf","CGf","CGr",
               "CCGG")
    # colls <- c(colls, "mq")
    df <- df_summary[, .SD, .SDcols=colls]

    # clean up
    rm(df_summary)
    colnames(df)[c(1)] <- c("CHR")
  
    # Remove loci with 0 distinguishable ref/alt alleles & no methylation
    df[, flag := ifelse(!is.na(ref)&total_counts>0, TRUE, ifelse(!is.na(total_counts_m)&total_counts_m>=3, TRUE,FALSE))]
      if(!is.na(table(df$flag)["FALSE"])) {
          warning(paste0("Removing loci with zero SNP allele counts and/or insufficient methylation - ", 
                         sum(table(df$flag)["FALSE"]), " positions removed.", sep = ""))
          df <- df[flag == TRUE,]
        }
    df[, flag := NULL]
    
    # Get Ensembl chrom annotation for compatibility with ASCAT
    if(UCSC == TRUE){df[,chrom := as.character(substr(df$CHR, 4, 5))]}
    if(UCSC == FALSE){ df[,chrom := CHR] ; df[,CHR := paste0("chr", as.character(CHR))]}
    df[,chrom := factor(chrom,levels=c(1:22,"X","Y"),ordered = TRUE)]
    
    # Compute methylation rate for CpGs
    df[, m := ifelse(width == 1, NA, ifelse(total_counts_m == 0, NA, M / (M + UM)))]
    # # print diganostics
    # cat("Bulk methylation rate calculation: complete!", sep = "\n", file = file_stdout_2, append = TRUE)
    
    # Compute BAF for SNPs
    df[, BAF := ifelse(is.na(ref), NA, ifelse(total_counts == 0, NA, alt_counts / total_counts))]
    # # print diganostics
    # cat("BAF calculation: complete!", sep = "\n", file = file_stdout_2, append = TRUE)
  
    # Flag tetranucleotide or dinucleotide positions which are duplicated due to SNP.
    df[, ID := paste(CHR,POS,sep="-")]
    dups <- df[!is.na(POS)&(duplicated(ID)|duplicated(ID,fromLast=TRUE)),]
    dups[, flag := ((BAF<=0.1)|(BAF>=0.9))]
    
    # If a position is duplicated due to a SNP overlapping with 2 CG/CCGG loci
    # check if the SNP is homozygous. If so, only one of the 2 loci can be a CG or CCGG 
    # Duplicate SNPs will not be removed at this stage at hetorozygous CGG/CCG
    if(dups[,sum(flag, na.rm=T)>0]){
      hom_dups <- dups[flag==TRUE,] 
      hom_dups[, flag := NULL]
      
      hom_dedup <- hom_dups %>% group_by(ID) %>% 
                   filter(total_counts_m == max(total_counts_m,na.rm=TRUE) & length(na.omit(M+UM))>0) %>%
                   arrange(factor(chrom,levels=c(1:22,"X","Y"),ordered = TRUE),start,end)
        
      df <- df[!(ID%in%hom_dups$ID),]
      df <- bind_rows(df, hom_dedup)
      df <- df[order(factor(chrom,levels=c(1:22,"X","Y"),ordered = TRUE),start,end),]
      rm(hom_dedup,hom_dups)
    };rm(dups)
    
    # format data.frame
    df[, ID := NULL]
    
    # Return empty BAM file if no data
    # If not, rownames() call below will raise error
    if(nrow(df) == 0){
      df_names = c("CHR","chrom","start","end","width",
                   "POS","ref","alt","alt_counts","ref_counts",
                   "total_counts","BAF","total_depth",
                   # "other_counts","all_counts",
                   "M","UM","total_counts_m","m",
                   # allele counts breakdown may be obtained by uncommenting this line + line 441.
                   #"Af","Ar","Cf","Cr","Tf","Tr","Gf","Gr","CAr","TGf","CGf","CGr", 
                   "CCGG")
      df <- data.frame(matrix(NA, nrow=1, ncol=length(df_names)))
      names(df) <- df_names
      df <- df[-1,] # Return empty dataframe
      return(df)
    } 
    
    # Format row names and col classes
    rownames(df) <- 1:nrow(df)
    if(getOption("scipen")==0){options(scipen = 999)}
    #df <- df[, lapply(.SD, as.character), by = 1:nrow(df)]
    
    # Ensure there are no NAs before saving
    flag <- df[,!is.na(CHR)]
    if(sum(flag==FALSE)>0){df <- df[flag,]};rm(flag)
  
    # Select columns for saving and put in a logical order
    df <- data.frame(df, stringsAsFactors = FALSE)
    df <- df[,c("CHR","chrom","start","end","width",
                "POS","ref","alt","alt_counts","ref_counts",
                "total_counts","BAF","total_depth",
                # "other_counts","all_counts",
                "M","UM","total_counts_m","m",
                # allele counts breakdown may be obtained by uncommenting this line + line 441.
                #"Af","Ar","Cf","Cr","Tf","Tr","Gf","Gr","CAr","TGf","CGf","CGr", 
                "CCGG")]
                # "mq")]
    
  return(df)
}
 
  # get_reads in parrallel
  if(n_cores>1){
    print(n_cores)
    
    # Set the cluster
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Pass libPaths to each worker on the cluster if incorrect libPaths in Rprofile.
       if(file.exists("~/.rprofile")|file.exists("~/.Rprofile")){
         x <- .libPaths()
         suppressWarnings(try(source("~/.rprofile"), silent=TRUE))
         suppressWarnings(try(source("~/.Rprofile"), silent=TRUE))
         y <- .libPaths()
         if(sum(!x%in%y)>0){
           clusterCall(cl, function(x) .libPaths(x), .libPaths())
         }
       }

    # run get reads function
    df_merged <- NULL
    df_merged <- foreach(j = 1:n_cores, .combine='rbind.data.frame', 
                         .packages=c("magrittr","Rsamtools", "GenomicRanges", "S4Vectors", 
                                     "IRanges", "dplyr", "data.table"),
                         .multicombine = TRUE) %dopar% {
                 df_Jth_core <- get_reads(i=i,j=j,n_cores=n_cores,bam_file=bam_file,
                                          segments_subset=segments_subset,mq=mq,test=test)
                                 }
    stopCluster(cl)
  }
  
  if(n_cores==1){
    df_merged <- get_reads(i=i,j=NULL,n_cores=NULL,bam_file=bam_file,
                          segments_subset=segments_subset,mq=mq,test=test)
    }
  invisible(gc())
  
  # Create file
  f_nm <- file.path(path_output, paste(patient_id, ".", sample_id, ".", i, ".SNPs.CpGs.fst", sep = ""))
  fst::write_fst(df_merged,  f_nm)
  cat(paste0("Written to: ", f_nm, "\n"))
}

# END
