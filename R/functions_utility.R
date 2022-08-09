#' Generate a filename that follows naming convention
#'
#' @param study.name character string; organizational name for a sequence of experiments
#' @param f.names character vector; file names/file paths that include study.name
#' @param prefix.name character string; unique identifier for the filename
#' @param suffix.name character string; unique identifier for the filename
#' @param f.extension character string; filename extension
#'
#' @return a filename that follows established naming convention; used to generate filenames for saving outputs
#' @export
#'
generate.out.name <- function(study.name=NULL,f.names=NULL,prefix.name=NULL,suffix.name=NULL,f.extension=".rds"){
  f.names.split <- strsplit(f.names,"/")
  exps.unique.names <- unique(sapply(f.names.split,function(i) grep(study.name,i,value = T)))
  exps.unique.names.split <- strsplit(exps.unique.names,"_")
  exps.unique.numbers <-  sapply(exps.unique.names.split,function(i){grep("^[0-9]{3}$",strsplit(i,"_"),value = T)})
  #
  out.name <- paste(study.name,
                    min(exps.unique.numbers),
                    "thru",
                    max(exps.unique.numbers),
                    sep = "_"
  )
  #
  if(!is.null(prefix.name)){
    out.name <- paste(prefix.name,out.name,sep="_")
  }
  if(!is.null(suffix.name)){
    out.name <- paste(out.name,suffix.name,sep="_")
  }
  #
  out.name <- paste(out.name,Sys.Date(),sep="_")
  out.name <- paste0(out.name,f.extension)
  return(out.name)
}
