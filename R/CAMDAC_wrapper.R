##  CAMDAC wrapper
##  Latest update: 04/05/2021
##  Version 1.0.0
##  Author: Elizabeth Larose Cadieux

# Load library
x <- c("optparse","stringr","GenomicRanges","Rsamtools","ASCAT","scales","ggplot2",
       "gridExtra","data.table","readr","fst","dplyr","parallel","doParallel")
# library(pryr) # only needed if mem_used() is used to determine object RAM.
invisible(lapply(x, function(y) suppressWarnings(suppressMessages(library(y, character.only = TRUE,
                    quietly=TRUE,warn.conflicts = FALSE)))))
rm(x)

# parse options
option_list = list(
  make_option(c("--patient_id"), type="character", default=NULL,
              help="Character string with patient id. [example : -patient_id = \'PX1\']"),
  make_option(c("--sample_id"), type="character", default=NULL,
              help="Character string with sample id. This must be filled in for both multi-region and single sample datasets. [example : -sample_id = \'T1\']"),
  make_option(c("--sex"), type="character", default="XY",
              help="Sex of the patient expressed as \'XX\' for female or \'XY\' for male [default  \"%default\"]."),
  make_option(c("--normal_infiltrates_proxy_id"), type="character", default=NULL,
              help="Input the patient-matched tumour-adjacent sample id as a string."),
  make_option(c("--normal_origin_proxy_id"), type="character", default=NULL,
              help="Input tissue-matched control sample id as a string for differential methylation analysis."),
  make_option(c("--patient_matched_normal_id"), type="character", default=NULL,
              help="Input patient-matched normal sample id as a string for genotyping purposes. \n
                    If no patient-matched normal sample is provided, heterozygous SNPs will be approximated from the tumour sample."),
  make_option(c("--bam"), type="character", default=NULL,
              help="Character string with full bam file name and path. \n
                    [example : -bam =\'/home/user/project_id/patient_id/bams/bam_file_sample1.bam\'] \n
                    Bam file must be sorted and indexed. \n
                    Index file must have the same file name and path with added .bai extension."),
  make_option(c("--ref"), type="character", default=NULL,
              help="Reference genome build. Either \'hg19\', \'GRCH37\',  \'hg38\' or  \'GRCH38\'. \n
                    The reference build MUST be one of these four options, but you may set build = NULL \n
                    if you are unsure which of these 4 builds your data is aligned to and let CAMDAC will detect the build for you."),
  make_option(c("--path_to_CAMDAC"), type="character", default=NULL,
              help="Character string containting the path to the CAMDAC library \n
                  [example : -ptc =\'/home/user/libs/CAMDAC/\']"),
  make_option(c("--wdir"), type="character", default="./",
              help="Path to to working directory [example : -wd =\'/home/user/project_id/\']"),
  make_option(c("--nc"), type="integer", default=1,
              help="Number of cores for parrallel processing [default %default].\n
                    We recommend running on 8-12 cores.")
); opt = parse_args(OptionParser(option_list=option_list))


# Store variables
patient_id = ifelse(is.null(opt$patient_id), error("You must input a patient id"), opt$patient_id)
sample_id = ifelse(is.null(opt$sample_id), error("You must input a sample id"), opt$sample_id)
sex = ifelse(is.null(opt$sex) , error(paste("You must input sex as \'XX\' for females or",
                                            "\'XY\' for males.", sep= " ")), opt$sex)
build = ifelse(is.null(opt$ref), error("You must input either \'hg19\', \'GRCH37\',  \'hg38\' or  \'GRCH38\'"),
               opt$ref)
path_to_CAMDAC=ifelse(is.null(opt$path_to_CAMDAC), error(paste("You must input the path to the CAMDAC library such as",
                                                               "\'/home/user/libraries/CAMDAC/\'", sep=" ")), opt$path_to_CAMDAC)
path = ifelse(is.null(opt$wdir),  warning(paste("You must input a path to the working directory where the tumour and",
                                                "control sample outputs will be stored such as",
                                                "\'/home/user/project_id/CAMDAC_results/\'", sep="\n")), opt$wdir)
# path/to/patient/dir (example: $HOME/project_id/)
n_cores = as.numeric(opt$nc)

# Process normal sample flags 
normal_infiltrates_proxy_id = ifelse(is.null(opt$normal_infiltrates_proxy_id), 
                                     error("You must provide a proxy for the normal tumour-infiltrating cells."), 
                                     opt$normal_infiltrates_proxy_id)
normal_origin_proxy_id = ifelse(is.null(opt$normal_origin_proxy_id) , 
                                error("You must provide a proxy for the normal cell of origin."), 
                                opt$normal_origin_proxy_id)
if(is.null(opt$patient_matched_normal_id)){
  cat("No patient-matched normal provided, heterozygous SNPs will be dertermined directly from the tumour sample.", 
  "LogR will be derived using a sex-matched panel of normal lung RRBS data from the TRACERx study.", 
  "Alternatively, users can build and provide their own sex-matched reference panel.", sep="\n")
  build_tmp = ifelse(build %in% c("hg19", "GRCH37"), "hg19", NULL)
  if(is.null(build_tmp)){
    error(paste0("RRBS coverage example reference panel only available for hg19 / GRCH37.", 
          " Please provide a normal sample for LogR calculation and copy number profiling.",
          " You may modify this wrapper to accept custom reference panels."))
  }
  reference_panel_coverage = file.path(path_to_CAMDAC, "pipeline_files", "example_normal_data", 
                                       paste0("normal_lung_cov_",build_tmp,"_", sex, ".fst"))
} else {
  patient_matched_normal_id = opt$patient_matched_normal_id
  reference_panel_coverage = NULL
}

# Set additional internal normal flags
is_normal = ifelse(sample_id %in% c(normal_infiltrates_proxy_id, normal_origin_proxy_id, patient_matched_normal_id), TRUE, FALSE)
is_patient_matched_normal = ifelse(sample_id == patient_matched_normal_id, TRUE, FALSE)

# Check that control sample has been run already
f_control1 <- file.path(path, patient_id, "Allelecounts", normal_infiltrates_proxy_id, paste(patient_id, normal_infiltrates_proxy_id,
                    "SNPs.CpGs.all.sorted.RData", sep="."))
f_control2 <- file.path(path, patient_id, "Allelecounts", normal_origin_proxy_id, paste(patient_id,  normal_origin_proxy_id,
                    "SNPs.CpGs.all.sorted.RData", sep="."))
if(is_normal == FALSE & (!file.exists(f_control1) | !file.exists(f_control2))) {
  stop(paste0("You must run CAMDAC on your normal sample(s) first.",
               "Do not change output directory structure and file names."))
}
rm(f_control1, f_control2)

# Store bam file name and carry out file checks
bam_file = ifelse(is.null(opt$bam), warning("You must input the sample bam file name and path"), opt$bam)

# Check that the bam_file exists
if(!file.exists(bam_file)) {
  stop(paste0("No bam file with name ", bam_file))
}

# Check that the index file is bam_file.bai
bai_file <- paste0(bam_file, ".bai")
if(!file.exists(bai_file)) {
  stop(paste0("No index file with name ", bai_file))
} ; rm(bai_file)

# Checkpoint
cat("Patient data output directory set to : ", path, "\n",
    "Do not alter directory structure or filenames", "\n", sep="")
  
# Source script to get SNP allele counts and CpG methylation data
source(paste0(path_to_CAMDAC,"/R/get_allele_counts.R"))
# Source script to format output and QC fragment size distribution
source(paste0(path_to_CAMDAC,"/R/format_output.R"))
# Source script to get copy number and purity estimates
source(paste0(path_to_CAMDAC,"/R/run_ASCAT.m.R"))
# Source script to QC and filter methylation data
source(paste0(path_to_CAMDAC,"/R/run_methylation_data_processing.R"))
# Source script to get purified tumour methylation rates, DMP and DMR calls
source(paste0(path_to_CAMDAC,"/R/get_pure_tumour_methylation.R"))
# Source script to get tumour-normal DMP and DMR calls
source(paste0(path_to_CAMDAC,"/R/get_differential_methylation.R"))

# Get SNP and CpG methylation (di)nucleotide counts
# Outputs (di)nucleotide counts, methylation rates at CpGs and BAF at SNPs.
for(a in 1:25){
    get_allele_counts(i=a, patient_id=patient_id, sample_id=sample_id, sex=sex, bam_file=bam_file,
                      mq=0, path=path, path_to_CAMDAC=path_to_CAMDAC, build=build, n_cores=n_cores, 
                      test=FALSE)
}

# Format output for copy numner and methylation analysis. Outputs a GRanges object.
format_output(patient_id=patient_id, sample_id=sample_id, sex=sex, is_normal=is_normal, 
              path=path, path_to_CAMDAC=path_to_CAMDAC, build=build, txt_output=FALSE)
                                           
# Get copy number (includes logR gc correction)
run_ASCAT.m(patient_id=patient_id, sample_id=sample_id, sex=sex,
            patient_matched_normal_id=patient_matched_normal_id,
            path=path, path_to_CAMDAC=path_to_CAMDAC, build=build, 
            min_normal=10, min_tumour=1,
            n_cores=n_cores, reference_panel_coverage=reference_panel_coverage)

# Process methylation info for copy number profiling and plot summary.
run_methylation_data_processing(patient_id=patient_id, sample_id=sample_id,
                                normal_infiltrates_proxy_id=normal_infiltrates_proxy_id, 
                                normal_origin_proxy_id=normal_origin_proxy_id,
                                path=path, min_normal=10, min_tumour=3, n_cores=n_cores, 
                                reference_panel_normal_infiltrates=NULL,
                                reference_panel_normal_origin=NULL)

if(!sample_id %in% c(normal_infiltrates_proxy_id, normal_origin_proxy_id)){
    # Get purified methylation rate 
    get_pure_tumour_methylation(patient_id=patient_id, sample_id=sample_id, sex=sex,
                                normal_infiltrates_proxy_id=normal_infiltrates_proxy_id, 
                                path=path, path_to_CAMDAC=path_to_CAMDAC, build=build,
                                n_cores, reseg=FALSE)

    # Get DMP and DMR calls
    get_differential_methylation(patient_id=patient_id, sample_id=sample_id,sex=sex,
                                 normal_origin_proxy_id=normal_origin_proxy_id,
                                 path=path,path_to_CAMDAC=path_to_CAMDAC,build=build,
                                 effect_size=0.2,prob=0.99,
                                 min_DMP_counts_in_DMR=5,min_consec_DMP_in_DMR=4,
                                 n_cores=n_cores, reseg=FALSE, bulk=FALSE)
}

# END
