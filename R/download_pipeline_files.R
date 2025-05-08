
#' Download CAMDAC pipeline files
#' 
#' @description CAMDAC pipeline files are required for analysis. This function downloads the files to
#' the output directory and unpacks them. By default, CAMDAC searches for the files in the
#' environment variable CAMDAC_PIPELINE_FILES. If this is missing, the current directory is used.
#' @param assay Sequencing assay. Either wgbs or rrbs.
#' @param directory Optional. Directory to download files to.
#' @export
download_pipeline_files <- function(bsseq, directory=NULL, quiet=TRUE){
  stopifnot(bsseq %in% c("wgbs", "rrbs", "test"))
  
  # Get download URL from CAMDAC index file
  url_index_file = system.file("extdata", "pipeline_files_urls.txt", package = "CAMDAC")
  urls = read.table(url_index_file, header=F, stringsAsFactors = F)
  names(urls) = c("bsseq", "link")
  link = urls[ urls$bsseq == bsseq, ][[2]]
  
  # Get download location
  #   If a directory is passed to the function, install there.
  if(!is.null(directory)){
    location=fs::path_expand(directory)
    fs::dir_create(location) # Ensure location exists
  } else {
  #   Else, get pipeline files location from environment variable
  #   The currect directory is used if environment variable is empty
    cpf_env = Sys.getenv("CAMDAC_PIPELINE_FILES")
    location = ifelse(cpf_env=="", ".", cpf_env)
    location = fs::path_expand(location)
  }
  
  # Ensure download directory path exists
  if(!fs::dir_exists(location)){
    fs::dir_create(location)
  }
  
  # Download pipeline files and unzip
  tf = tempfile()
  tryCatch({
    download.file(link, destfile=tf, method="wget", quiet=quiet)
  }, error=function(e){
    logerror("Pipeline files for {bsseq} could not be downloaded from {link}.")
    stop()
  })
  untar(tf, exdir=location)
  
  loginfo("Pipeline files for {bsseq} downloaded to {location}")
  return(location)
}
