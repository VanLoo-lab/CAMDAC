.onLoad <- function(libname, pkgname){
  # Silence load warnings for conflicting global methods from imports.R packages
  options(conflicts.policy = list(
    error=FALSE,
    warn=FALSE,
    generics.ok=TRUE,
    depends.ok=TRUE,
    can.mask=c("select")
  ))

  # Setup CAMDAC package logger
  library(logging)
  logging::logReset() # Fixes root logger, but does this mess up other packages?
  camdac_logger <- logging::getLogger(pkgname)
  logging::setLevel("INFO", camdac_logger)
  logging::addHandler(logging::writeToConsole,
                      logger=pkgname,
                      formatter = function(record) {
                        sprintf("[%s] | %s | %s | %s",
                                record$timestamp,
                                record$logger,
                                record$levelname,
                                record$msg)
                      })
  logging::loginfo(sprintf("Package loaded."), logger="CAMDAC")
}