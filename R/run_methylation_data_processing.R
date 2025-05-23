##  Process raw CpG methylation info
##  Latest update: 04/05/2021
##  Version 1.0.0
##  Author: Elizabeth Larose Cadieux

#' Filter bulk tumour and normal methylation data, get methylation rate highest density interval (HDI)
#' and plot raw methylation info
#' \code{run_methylation_data_processing}
#'
#' @param patient_id Character variable containting the patient ID
#' @param sample_id Character variable with the (control or tumour) sample ID
#' @param normal_infiltrates_proxy_id Character variable with the sample ID of
#' the tissue-matched normal acting as proxy for the tumour infiltrating
#' normal cells. Ideally, this is a patient and tissue-matched tumour adjacent normal sample.
#' @param normal_origin_proxy_id Character variable with the sample ID
#' of the normal to be used as a proxy for the tumour cell of origin in
#' differential methylation analyses.
#' @param path Character path variable pointing to the desired working directory.
#' This is where the output will be stored.
#' @param min_tumour Numerical value correspdonding to the minimum counts threshold
#' in the tumour sample CpGs inclusion
#' @param min_normal Numerical value correspdonding to the minimum counts threshold for
#' the normal CpGs to be included
#' @param n_cores Numerical value correspdonding to the number of cores for parallel processing
#' @param reference_panel_normal_infiltrates Default is NULL. Character string with the complete
#' path to a reference methylation profile for the tumour normal infiltrates as a .fst file.
#' @param reference_panel_normal_origin Default is NULL. Character string with the complete
#' path to your reference methylation profile for the tumour cell of origin as a .fst file.
#' 
#' If a patient-matched proxy for the normal infiltrates and/or the normal cell of origin is not
#' available, a reference panel may be constructed from different individuals and used as a substitute.
#'
#' The reference samples should be at the very least sex-matched.
#'
#' The reference should be saved as a .fst file with the following columns:
#' CHR  start      end        M_n        UM_n       m_n        cov_n
#' <factor>    <integer>  <integer>  <numeric>  <numeric>  <numeric>  <numeric>
#'
#' where each row is a CpG or CCpGG with coordinates CHR:start-end 
#' The start and end columns correspond to the 5'-C and 3'-G coordinate, respectively.
#' M_n is the number of reads supporting of the methylated allele
#' UM_n is the number of reads supporting of the unmethylated allele
#' m_n is the normal methylation rate (M_n / (M_n+UM_n))
#' cov_n is the total CpG methylation informative reads counts (M_n+UM_n)
#'
#' @return GRanges object in .RData file
run_methylation_data_processing <- function (patient_id,sample_id,
                                             normal_infiltrates_proxy_id,
                                             normal_origin_proxy_id,
                                             path,min_normal=10,min_tumour=3,n_cores,
                                             reference_panel_normal_infiltrates=NULL,
                                             reference_panel_normal_origin=NULL){
    
if(getOption("scipen")==0){options(scipen = 999)}
# define normals
n_i_id <- normal_infiltrates_proxy_id
n_o_id <- normal_origin_proxy_id
normal_ids <- unique(c(n_i_id,n_o_id))

# define column name suffix for normal and bulk tumour samples
# this enables to distinguish between:
# the proxy supplied for the normal infiltrates (_n_i) for use in deconvolution
# the normal cell of origin (_n_o) for use in DMP/DMR calling
# and the bulk tumour (_b)
if(sample_id %in% normal_ids){
    suffix <- ifelse(n_i_id==n_o_id, "_n", ifelse(sample_id==n_i_id, "_n_i", "_n_o"))
} else {
    suffix <- "_b"
}

# Set paths
path_patient <- file.path(path, patient_id)
dir.create(file.path(path_patient, "Methylation",sample_id), recursive=TRUE,showWarnings=FALSE)
path_output <-paste0(path_patient, "/Methylation/",sample_id,"/")

# Load data
f_nm = paste0(path_patient,"/Allelecounts/", sample_id, "/",
              patient_id, ".", sample_id, ".SNPs.CpGs.all.sorted.RData")
load(f_nm)
rm(f_nm)

# Identify CpGs that meet the minimum coverage treshold
if(sample_id %in% normal_ids){dt <- dt_normal ; rm(dt_normal) ; min<-min_normal}
if(!sample_id %in% normal_ids){dt <- dt_tumour ; rm(dt_tumour) ; min<-min_tumour}

# Get non-numereric column IDs
k1 <- which(colnames(dt)=="alt_counts") 
k2 <- which(colnames(dt)=="CCGG") 
vec <- k1:k2
flag <- dt[1, lapply(.SD, class), .SDcols=c(vec)]%in%c("numeric", "integer")

# If applicable, convert relevant columns to numerical
if(length(vec[!flag])>0){
dt[, as.character(colnames(dt)[vec[!flag]]) := lapply(.SD, as.numeric), .SDcols=vec[!flag]]
}

# Select relevant columns
dt1 <- dt[, .(CHR, start, end, POS, ref, alt)]
dt2 <- dt[, .SD, .SDcols=vec]
dt <- cbind(dt1, dt2)
rm(dt1, dt2, k1, k2, vec, flag)

# filter out SNPs outside CpG / CCGG and CpG with low coverage
flag <- dt[, (end-start) > 0 & !is.na(m) & total_counts_m >= min]
dt[flag, mloci := 1:sum(flag)]
dt_sample_m <- dt[!is.na(mloci),]
rm(min,flag)

if(!sample_id %in% normal_ids){
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
dt_sample_m[, offset := ifelse((end-start)>2,1,0)] 
dt_sample_m[, ID := 1:nrow(dt_sample_m)]
tmp <- dt_sample_m[, .(ID, CHR, start=start+offset, end=end-offset)]
setkeyv(tmp, cols=c("CHR", "start", "end"))
ov <- foverlaps(dt_SNPs, tmp, nomatch = 0)

# store label and tumour BAF for downstream purified tumour methylation rate calculation
dt_sample_m[, c("BAF") := ov[match(dt_sample_m$ID,ov$ID), .(BAF)]]
dt_sample_m[, c("offset", "ID") := NULL]  
rm(dt_SNPs, tmp)
}

# Assign CpG class
dt_sample_m[, class := factor(ifelse(!is.na(BAF) & CCGG==0 & !is.na(BAF),"SNP CpG",
                              ifelse(!is.na(BAF) & CCGG>0 & !is.na(BAF) &
                                     POS%in%c(start+1, end-1),"SNP CCGG",
                              ifelse(CCGG==0, "non-SNP CpG",
                              ifelse(CCGG>0, "non-SNP CCGG", NA)))),
                        levels = c("SNP CpG", "SNP CCGG", "non-SNP CpG", "non-SNP CCGG"), 
                        ordered = TRUE)]

# Make plots
df_sample_m <- data.frame(dt_sample_m, stringsAsFactors = F)
outfile <- paste0(path_output,patient_id,"_", sample_id,"_methylation_rate_summary.pdf")
plot_methylation_info(df=df_sample_m,outfile=outfile)
rm(df_sample_m)

# Checkpoint
if(!sample_id%in%normal_ids){cat("Bulk methylation rates plots generated\n")}
if(sample_id%in%normal_ids){cat("Normal methylation rates plots generated\n")}

# Format methylation into Granges and apply mim CpG coverage treshold for DMP calculation
dt_sample_m <- format_methylation_df(dt=dt_sample_m,sample_id=sample_id,normal_id=normal_ids,
                                     suffix=suffix,path_output=path_output,n_cores=n_cores)

# Assign matched normal and compute HDI
if(!sample_id%in%normal_ids){
  # set control filename(s) and reference panels if relevant
  filename_control <- paste0(path_patient, "/Methylation/",unique(normal_ids),"/dt_normal_m.RData")
  reference_panel <- c(reference_panel_normal_infiltrates, reference_panel_normal_origin)
  if(sum(!is.null(reference_panel))==0){
    reference_panel <- NULL
  }
  if(length(filename_control)>1 & sum(!is.null(reference_panel))>0){
      reference_panel <- character(length=2)
      reference_panel[1] <- ifelse(!is.null(reference_panel_normal_infiltrates), reference_panel_normal_infiltrates, NA)
      reference_panel[2] <- ifelse(!is.null(reference_panel_normal_origin), reference_panel_normal_origin, NA)
  }
 
  # set tumour data.table key
  setkeyv(dt_sample_m,  cols = c("CHR", "start", "end"))
  dt_normal_m <- NULL

  # Add matched normal info
  for(i in 1:length(filename_control)){
     if(file.exists(filename_control[i])){ 
      load(filename_control[i])
     } else if(sum(!is.null(reference_panel))>0 & !is.na(reference_panel[i])){
      # If you don't have patient matched normal
      dt_normal_m <- read_fst(path=reference_panel[i], as.data.table=TRUE)

      # Ensure all required columns are present
      if(!is.null(dt_normal_m)&sum(c("CHR", "start", "end", "M_n", "UM_n", "m_n", "cov_n")%in%colnames(dt_normal_m))!=7){
         stop(paste0("Required columns not all present in reference panel for the ",
                     c("tumour-infiltrating normal cells", "cell of origin")[i]))
      }
      # Ensure there are no NAs
      dt_normal_m <- dt_normal_m[!is.na(m_n),]

      # Calculate HDI for the methylation rate
      n<- nrow(dt_normal_m)
      M=dt_normal_m$M_n;UM=dt_normal_m$UM_n
      vec <- cbind.data.frame(low=numeric(length=n), high=numeric(length=n))
      vec[,1:2]<- do.call(rbind, parallel::mclapply(1:n, function(i,M,UM) 
          HDIofICDF(ICDFname=qbeta, credMass=.99, shape1=M[i]+1, shape2=UM[i]+1),
                                        mc.cores=n_cores, M=M, UM=UM))
      # Checkpoint
      cat("Reference profile methylation rates 99% highest density intervals annotated\n")

      # Add HDI to data.table
      dt_normal_m[, as.character(c("m_n_low", "m_n_high")) := as.list(vec)]
      rm(n,M,UM,vec)

      # re-order and rename columns
      dt_normal_m <- dt_normal_m[, .SD, .SDcols=c("CHR", "start", "end", "M_n", "UM_n", "m_n", "cov_n", "m_n_low", "m_n_high")]
      if(length(filename_control)>1&i==1){
         colnames(dt_normal_m)[4:ncol(dt_normal_m)] <- 
            c(paste0(c("M_", "UM_", "m_", "cov_"), "n_i"), "m_n_i_low", "m_n_i_high")
      } 
       if(length(filename_control)>1&i==2){
          colnames(dt_normal_m)[4:ncol(dt_normal_m)] <- 
            c(paste0(c("M_", "UM_", "m_", "cov_"), "n_o"), "m_n_o_low", "m_n_o_high")
      } 
  
    } else {
      stop(paste0("No proxy provided for the normal ", c("tumour- infiltrating cells", " cell of origin")[i]), ".\n",
                  "You must provide a proxy normal profile for the tumour infiltrating cells for tumour deconvolution and\n
                  a proxy for the cell of origin for differential methylation analysis calling.")
    }

    # Check seqlevels are matched
    if(!grepl("chr", dt_sample_m$CHR[1])){
        dt_sample_m[, CHR := paste0("chr", CHR)]
    }
    if(!grepl("chr", dt_normal_m$CHR[1])){
        dt_normal_m[, CHR := paste0("chr", CHR)]
    }

    # set normal key
    setkeyv(dt_normal_m,  cols = c("CHR", "start", "end"))
 
    # Annotate normal data on tumour data.table
    overlaps <- foverlaps(dt_normal_m, dt_sample_m, nomatch=0, type="equal")
    if(length(filename_control)==2 & i==1){
       coll <- colnames(dt_sample_m)
       coll_n <- colnames(dt_normal_m)[!grepl("(CHR)|(start)|(end)", colnames(dt_normal_m))]
       dt_sample_m <- overlaps[, .SD , .SDcols = c(coll, coll_n)]
       
       # set key
       setkeyv(dt_sample_m,  cols = c("CHR", "start", "end"))
     }
    if(length(filename_control)==2 & i==2){
       coll_n2 <- colnames(dt_normal_m)[!grepl("(CHR)|(start)|(end)", colnames(dt_normal_m))]
       dt_tumour_and_normal_m <- overlaps[, .SD, .SDcols = c(coll, coll_n, coll_n2)]
    }
    if(length(filename_control)==1){
       coll <- colnames(dt_sample_m)
       coll_n <- colnames(dt_normal_m)[!grepl("(CHR)|(start)|(end)", colnames(dt_normal_m))]
       dt_tumour_and_normal_m <- overlaps[, .SD , .SDcols = c(coll, coll_n)]
     }
  }
  rm(overlaps, dt_normal_m)
  
  # to fix later
  dt_tumour_and_normal_m <- unique(dt_tumour_and_normal_m)
  
  # Save
  save(dt_tumour_and_normal_m, file=paste0(path_output, "dt_tumour_and_normal_m.RData"))
  }
}

#' Format methylation rates
#' \code{format_methylation_df}
#'
#' @param dt data.table containing the methylation information for each CpG
#' @param sample_id sample ID
#' @param normal_ids sample ID of normal sample(s)
#' @param path_output output directory
#' @param n_cores number of threads for HDI calculation
#' @param suffix string containing the column names suffix for normal samples
#'        This is to distinguish between the proxy supplied for the normal infiltrates
#'        for use in deconvolution and the normal cell of origin for use in DMP/DMR calling
#' @param trim Logical value establishing whether regions with extremely high coverage be trimmed or not
#'
#' @return A GRanges object with all the CpG loci, their coverage, counts methylated and methylation rate
format_methylation_df <- function (dt,sample_id,normal_ids,path_output,n_cores,suffix,trim=FALSE) {
  
  # Get total cov (UM includes hetorozygous SNP non-CpG allele counts)
  dt[, cov := UM + M]
  cols <- c("CHR", "start", "end", "cov", "M", "UM", "m",
            "BAF", "POS", "ref", "alt")
  dt <- dt[, .SD, .SDcols=cols]
  rm(cols)
                       
  # Annotate heterozygous SNPs
  dt[, SNP := ifelse(!is.na(BAF) & BAF >= 0.15 & BAF <= 0.85, 1, 0)]
  # dt[, table(SNP)]

  # Remove positions tetranucleotide or dinucleotide which are duplicated 
  # This is due to 2 neighbouring SNP co-locolising with the CG/CCGG motif
  dt[, dups := length(cov)-1, by=.(CHR,start,end)]
  dups <- dt[dups>0,]
  
  # table(dups[,.(test=length(c(unique(M), unique(UM)))==2), 
  #               keyby = .(CHR,start,end)]$test)
  dedups <- dups[,.("M"=ifelse(sum(SNP==0)>0, # if there is a CG-destroying SNP
                       M[SNP==0][order(cov[SNP==0],decreasing=TRUE)][1],
                       M[order(cov,decreasing=TRUE)][1]),
                   "UM"=ifelse(sum(SNP==0)>0, # if there is a CG-destroying SNP
                       UM[SNP==0][order(cov[SNP==0],decreasing=TRUE)][1],
                       UM[order(cov,decreasing=TRUE)][1])), 
                  keyby = .(CHR,start,end)]
  dedups[, c("cov", "m") := list(M+UM, M/(M+UM))]
  
  # join deduplicated data
  dt <- dt[dups==0,] 
  dt[, c("dups", "BAF", "POS", "ref", "alt", "SNP") := NULL]
  dt <- rbind(dt, dedups)
  rm(dups,dedups)

  # order deduplicated data
  chr_names <- paste0("chr", c(1:22, "X", "Y"))
  dt <- dt[order(factor(CHR, levels=chr_names, ordered=TRUE), start, end)]
    
  # Order cols
  dt <- dt[, .SD, .SDcols=c("CHR", "start", "end", "M","UM","m","cov")] 
  dt <- dt[!is.na(m),] # ensure no NAs are introduce by dedup
  
  # Calculate HDI for the methylation rate
  n<- nrow(dt)
  M=dt$M;UM=dt$UM
  vec <- cbind.data.frame(low=numeric(length=n), high=numeric(length=n))
  vec[,1:2]<- do.call(rbind, parallel::mclapply(1:n, function(i,M,UM) 
                      HDIofICDF(ICDFname=qbeta, credMass=.99, shape1=M[i]+1, shape2=UM[i]+1),
                                mc.cores=n_cores, M=M, UM=UM))

  # Add HDI to data.table
  dt[, as.character(c("m_x_low", "m_x_high")) := as.list(vec)]
  rm(n,M,UM,vec)
 
  # format column names
  colnames(dt)[4:ncol(dt)] <-
    c(paste(c("M","UM","m","cov"), suffix, sep=""),
      paste0("m", suffix, "_low"),paste0("m", suffix, "_high"))

  if(sample_id%in%normal_ids){
    # Checkpoint
    cat("Normal methylation rates 99% highest density intervals annotated\n")

    # if required remove trim high coverage sites (probs = poor alignment)
    if(trim == TRUE){
      q <- dt[, quantile(cov, 0.9999)]
      dt <- dt[cov < q, ]
      rm(q)
    }
     
    # save dt
    dt_normal_m <- dt
    save(dt_normal_m, file=paste0(path_output, "dt_normal_m.RData"))
  }
  
  if(!sample_id%in%normal_ids){
    # Checkpoint
    cat("Bulk methylation rates 99% highest density intervals annotated\n")

    ## save dt
    #dt_tumour_m <- dt
    #save(dt_tumour_m, file=paste0(path_output, "dt_tumour_m.RData"))
  }
  
  return(dt)
}

# Arguments:
#'  @param ICDFname is R's name for the inverse cumulative density function
#' of the distribution.
#' @param credMass is the desired mass of the HDI region.
#' @param tol is passed to R's optimize function, 
#' the lower the tolerance,the longer the optimisation, but the higher the accuracy.
#' tol=1e-4 gives values of the same accurary as our max resolution
#' Return value:
#' Highest density iterval (HDI) limits in a vector.
#' Example of use: For determining HDI of a beta(30,12) distribution, type
#' HDIofICDF( qbeta , shape1 = 30+1 , shape2 = 12+1 )
#' Notice that the parameters of the ICDFname must be explicitly named;
#' e.g., HDIofICDF( qbeta , 30+1 , 12+1 ) does not work.
#' Adapted and corrected from Greg Snow's TeachingDemos package.

# Source fct outside of loop to speed up code
intervalWidth =  function(lowTailPr,ICDFname,credMass, ... ) {
ICDFname(credMass+lowTailPr, ... ) - ICDFname(lowTailPr, ... )
}

HDIofICDF = function(ICDFname, credMass=0.99 , tol=1e-4, ... ) {
  
  incredMass = 1.0 - credMass
  
  optInfo = optimize(f = intervalWidth, interval = c(0,incredMass) , ICDFname=ICDFname , credMass=credMass , tol=tol , ... )
  
  HDIlowTailPr = optInfo$minimum
  vec <- setNames(object = ICDFname(c(HDIlowTailPr, credMass+HDIlowTailPr), ... ), nm = c("low", "high"))
  return(vec)
}

#' Plot Methylation 
#' 
#' \code{plot_methylation_info} returns the df_sample with annotated q-value for each CpG
#'
#' @param df_sample data.frame with methylation info
#' @param outfile character srting with output pdf filename
#' 
#' @return pdf w/ methylation rate distribution, biases at polymorphic and non-polymorphic CG/CCGG and coverage distribution 
plot_methylation_info <- function (df_sample, outfile) {
  
  alph <- ifelse(df_sample$class %in% c("SNP CpG", "SNP CCGG"), "SNP", "non-SNP")
  p1 <- ggplot(data=df_sample, aes(x=m, y=after_stat(density), color=class, fill = class, alpha=alph)) + 
        ylab("Normalised density") + theme_classic() +
        geom_density(bw= 0.025) + 
        scale_x_continuous(name="CpG methylation rate", breaks=seq(0,1,.1)) +
        scale_color_manual(name = "", values = c("SNP CpG"="darkred","SNP CCGG"="darkblue", 
                           "non-SNP CpG"="lightsalmon","non-SNP CCGG"="lightblue")) +
        scale_fill_manual(name = "", values = c("SNP CpG"="darkred","SNP CCGG"="darkblue", 
                          "non-SNP CpG"="lightsalmon","non-SNP CCGG"="lightblue")) +
        scale_alpha_manual(name = "", values = c("SNP"=0.1,"non-SNP"=0.5)) +
        theme(legend.position="none") + ggtitle("A.") #+
 
  p2 <- ggplot(data=df_sample, aes(x=m, y=ggplot2::after_stat(count),color=class, fill = class)) + 
        theme_classic() + 
        geom_histogram(binwidth=0.025,alpha = 0.25) + 
        scale_x_continuous(name="CpG methylation rate", breaks=seq(0,1,.1)) +
        scale_y_continuous(name= "CpG counts",
                           breaks=seq(0,9,1)*10^5, 
                           labels=parse(text= paste0(seq(0,9,1),"%*%", 
                           scales::math_format(expr=10^.x)(c(5))))) +
        scale_fill_manual(name = "", values = c("SNP CpG"="darkred","SNP CCGG"="darkblue", 
                          "non-SNP CpG"="lightsalmon","non-SNP CCGG"="lightblue")) +
        scale_color_manual(name = "", values = c("SNP CpG"="darkred","SNP CCGG"="darkblue", 
                           "non-SNP CpG"="lightsalmon","non-SNP CCGG"="lightblue")) +
        theme(legend.position="none") + ggtitle("B.") #+theme(legend.text=element_text(size=8)) 
   
  tmp <- data.frame(table(df_sample$class))
  tmp <- tmp[order(tmp$Var1,decreasing = T),] ; colnames(tmp)[1:2] <- c("class", "frequency")
  tmp$anno_start <- c(0, cumsum(tmp$frequency[1:3])) ; tmp$anno_end <- cumsum(tmp$frequency)
  tmp$mid_point <- (tmp$anno_end-tmp$anno_start) / 2 + tmp$anno_start
  tmp$labels <- paste0(round(tmp$frequency/sum(tmp$frequency)*100, digits=0),"%")
  
  myTableGrob <- function(dt, title_v, fontsize_v = 10){
    #' Create custom table grob with title
    #' @description Creates table grob in format that is most common for my usage.
    #' @param dt Data.table that the grob will be made out of
    #' @param title_v Title for display
    #' @param fontsize_v Fontsize for title. Default is 14 (goes well with my_theme)
    #' 
    
    ## Table
    table_grob <- gridExtra::tableGrob(dt, rows = rep('', nrow(dt)), theme = ttheme_minimal(base_size=8,vjust=0, hjust=0))
    ## Title
    title_grob <- grid::textGrob(title_v, gp = grid::gpar(fontsize = fontsize_v),x=0,hjust=0)
    ## Add title
    table_grob <- gtable::gtable_add_rows(table_grob, heights = grid::grobHeight(title_grob) + unit(5,'mm'), pos = 0)
    table_grob <- gtable::gtable_add_grob(table_grob, title_grob, 1, 1, 1, ncol(table_grob), clip = "off")
  }
  
  df_sample_tmp <- data.table(df_sample)
  x <- df_sample_tmp[,.(y = median(total_depth)), by=.(class)]
  tmp$"median depth" <- x[match(tmp$class, x$class), ifelse(y>=0.5, ceiling(y), floor(y))]
  
  x <- df_sample_tmp[,.(y = median(m)), by=.(class)]
  tmp$"median m" <- x[match(tmp$class, x$class), 
      ifelse(as.integer(abs(y)*100) %% 10 >= 5, ceiling(y*100)/100, floor(y*100)/100)]
  rm(df_sample_tmp, x)
  
  p3 <- myTableGrob(tmp[,c("class", "frequency", "median depth", "median m")], "C.",fontsize=14)
  p4 <- ggplot(tmp, aes(x=factor(1), y=frequency, fill=class, color=class, label=labels)) + 
        theme_classic() + geom_bar(width = 1, alpha=0.25,stat = "identity") +  
        ggtitle("D.") + geom_text(size = 3, position = position_stack(vjust = 0.5), color="black") +
        theme(axis.text = element_blank(), axis.title = element_blank(), axis.line = element_blank(),
              axis.ticks = element_blank(), title = element_text(vjust=0,hjust=0)) +
        scale_fill_manual(name = "", values = c("SNP CpG"="darkred", "SNP CCGG"="darkblue", 
                                                "non-SNP CpG"="lightsalmon", "non-SNP CCGG"="lightblue")) +
        scale_color_manual(name = "", values = c("SNP CpG"="darkred", "SNP CCGG"="darkblue", 
                                                 "non-SNP CpG"="lightsalmon", "non-SNP CCGG"="lightblue")) + 
        coord_polar("y") 
  
# box_cov1 <- ggplot(data=df_sample, aes(x=class, y=total_depth, color = class, fill=class)) +
#   theme_classic() + xlab("CpG context") +
#   geom_violin(position=position_dodge(1)) + 
#   scale_y_log10(name= "log(Informative CpG coverage)",
#                 breaks = scales::trans_breaks("log10", function(x) 10 ^ x, n = 4)(c(1e1, 1e4)),
#                 labels= parse(text= unlist(scales::trans_format("log10", 
#                         scales::math_format(expr=10^.x))(c(1e1, 1e2, 1e3, 1e4))))) +   
#   scale_x_discrete(name = "CpG context", breaks = c("SNP CpG","SNP CCGG", "non-SNP CpG","non-SNP CCGG"), 
#                    labels = c( "SNP\nCpG", "SNP\nCCGG", "non-SNP\nCpG", "non-SNP\nCCGG"))  +
#   scale_color_manual(name = "", values = c("SNP CpG"="darkred", "SNP CCGG"="darkblue", 
#                                            "non-SNP CpG"="lightsalmon", "non-SNP CCGG"="lightblue")) +
#   scale_fill_manual(name = "", values = c("SNP CpG"="darkred", "SNP CCGG"="darkblue", 
#                                           "non-SNP CpG"="lightsalmon", "non-SNP CCGG"="lightblue")) +
#   theme(axis.ticks.x = element_blank(),axis.text.x=element_text(size=8)) +
#   ggtitle("E.") #+ theme(legend.position="none")
  p5 <- ggplot(data=df_sample, aes(x=total_depth, y=ggplot2::after_stat(count))) +
        theme_classic() + 
        geom_histogram(bins=50,alpha = 0.25, col="grey15") + 
        scale_x_continuous(name= "CpG coverage", 
                           breaks=seq(0,200,25), limits=c(1,200)) +
#        scale_x_log10(name= "log(Informative CpG coverage)",
#                      breaks = scales::trans_breaks("log10", function(x) 10 ^ x, n = 4)(c(1e1, 1e4)),
#                      labels= parse(text= unlist(scales::trans_format("log10", 
#                      scales::math_format(expr=10^.x))(c(1e1, 1e2, 1e3, 1e4))))) + 
         scale_y_continuous(name= "CpG counts",
                            breaks=seq(0,5,1)*10^5, 
                            labels=parse(text= paste0(seq(0,5,1),"%*%", 
                            scales::math_format(expr=10^.x)(c(5))))) +
         theme(axis.text.x=element_text(size=8)) + ggtitle("E.")

  Sys.sleep(1)
  
  gglist1 <- list(p1, p2)
  gglist2 <- list(p3, p4, p5)
  my_layout <- rbind(c(NA,rep(1,12),rep(3,12)),
                     c(NA,rep(1,12),rep(3,12)),
                     c(NA,rep(1,12),NA,rep(4,9), NA, NA),
                     c(NA,rep(1,12),NA,rep(4,9),NA, NA),
                     c(rep(2,13),NA,rep(4,9),NA, NA),
                     c(rep(2,13),rep(5,12)),
                     c(rep(2,13),rep(5,12)),
                     c(rep(2,13),rep(5,12)))
  Sys.sleep(1)
  gd <- arrangeGrob(grobs = c(gglist1,gglist2), layout_matrix = my_layout)
  Sys.sleep(1)
  
  ggsave(outfile,device = "pdf", gd, width=9, height = 6, units = "in")
  
  #clean
  suppressWarnings(rm(p1,box_cov1,p2,p3, p4, gd, my_layout, gglist1,gglist2, tmp))
}

# END
