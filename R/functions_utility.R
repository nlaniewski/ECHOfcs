#' Generate a filename that follows naming convention
#'
#' @param study.name character string; organizational name for a sequence of experiments
#' @param f.names character vector; file names that include study.name
#' @param prefix.name character string; unique identifier for the filename; what is it?
#' @param f.extension character string; filename extension
#'
#' @return a filename that follows established naming convention; used to generate filenames for saving outputs
#' @export
#'
generate.out.name <- function(study.name=NULL,f.names=NULL,prefix.name=NULL,f.extension=".rds"){
  f.names.split <- strsplit(f.names,"/")
  exps.unique.names <- unique(sapply(f.names.split,function(i) grep(study.name,i,value = T)))
  exps.unique.numbers <- stringr::str_extract(exps.unique.names,"[0-9]{3}")
  #
  paste0(paste(prefix.name,
               study.name,
               min(exps.unique.numbers),
               "thru",
               max(exps.unique.numbers),
               Sys.Date(),
               sep = "_"),
         f.extension)
}
