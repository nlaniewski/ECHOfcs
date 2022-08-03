#' Get .fcs marker names from .fcs file path(s)
#' @description Returns .fcs marker names based on established naming convention
#'
#' @param fcs.path .fcs file path
#' @param name.split character value used to split .fcs stain name: '$P##S'
#' @param split.position position of the stain name split that contains the marker name
#' @param selected.markers optional character vector of selected marker names; for getting '$P##N' column names
#'
#' @return character vector of marker names; or the column names of selected.markers
#' @export
#'
fcs.markers <- function(fcs.path,name.split="_",split.position=2,selected.markers=NULL){
  tmp.header <- flowCore::read.FCSheader(fcs.path)[[1]]
  ##
  p.names <- sapply(c("N","S"),function(i){
    p <- tmp.header[grep(paste0("P[0-9]+",i), names(tmp.header), value = T)]
    p <- p[order(as.numeric(stringr::str_extract(names(p),"[0-9]+")))]
  })
  ##
  ps.selected <- sapply(strsplit(p.names$S,name.split),'[',split.position)
  ps.selected <- ps.selected[!is.na(ps.selected)]
  ##
  if(!is.null(selected.markers)){
    selected.names <- sub("S","N",names(ps.selected)[which(ps.selected %in% selected.markers)])
    selected.markers.names <- unname(p.names$N[which(names(p.names$N) %in% selected.names)])
    return(selected.markers.names)
  }else{
    return(ps.selected)
  }
}

#' Check if markers exist in all .fcs files
#'
#' @param fcs.paths .fcs file paths
#' @param markers.vec a character vector of marker names
#'
#' @return logical TRUE or FALSE; message or stop
#' @export
#'
#'
fcs.markers.check <- function(fcs.paths,markers.vec,...){
  fcs.markers.list <- unique(lapply(fcs.paths,fcs.markers,...))
  if(Reduce(all,lapply(fcs.markers.list,function(i)all(markers.vec %in% i)))){
    message("All named markers found in .fcs files")
    return(TRUE)
  }else{
    stop("Marker conflict")
  }
}

#' Read .fcs file; keep only selected markers
#'
#' @param fcs.path .fcs file path
#' @param selected.markers character vector of marker names to keep
#' @param include.scatter logical; for fluor/flow based assays, include light-scatter channels
#'
#' @return a flowFrame; no transformation applied; no truncation
#' @export
#'
read.fcs.selected.markers <- function(fcs.path, selected.markers,include.scatter=F){
  selected.markers.names <- fcs.markers(fcs.path, selected.markers=selected.markers)
  if(include.scatter){
    selected.markers.names <- c(selected.markers.names,
                                paste(rep(c('FSC','SSC'),each=3),c('A','H','W'),sep = "-")
    )
  }
  fcs <- flowCore::read.FCS(fcs.path,
                            column.pattern = paste0(selected.markers.names,collapse = "|"),
                            transformation = F, truncate_max_range = F)
  ##
  return(fcs)
}

#' Read .fcs files; keep only selected markers; parallelized
#'
#' @param fcs.paths .fcs file paths
#' @param selected.markers character vector of marker names to keep
#'
#' @return a list of individual flowFrames
#' @export
#'
read.fcs.selected.markers.parallel <- function(fcs.paths,selected.markers){
  if(fcs.markers.check(fcs.paths,selected.markers)){
    message("Making parallel clusters; reading .fcs files")
    cl <- parallel::makeCluster(parallel::detectCores()-1)
    fcs.list <- parallel::parSapply(cl, fcs.paths, function(i,m.vec=selected.markers){
      fcs <- ECHOfcs::read.fcs.selected.markers(fcs.path=i,m.vec)
    })
    parallel::stopCluster(cl)
    #rm(cl)
    return(fcs.list)
  }
}

#' Generate some meta-data from .fcs file path/names; requires an established naming convention
#'
#' @param fcs.paths .fcs file paths
#' @param name.split character value used to split .fcs file name
#'
#' @return data.frame containing some meta-data; factored columns
#' @export
#'
mdat.frame.from.paths.echo <- function(fcs.paths,name.split="_"){
  if(class(fcs.paths)!="character"&!is.vector(fcs.paths)){
    stop("fcs.paths needs to be a character vector of .fcs file paths")
  }
  split.names.vec <- lapply(fcs.paths, function(i){
    batch <- sub("ECHO_","",stringr::str_extract(i,"ECHO_[0-9]{3}"))
    batch.date <- stringr::str_extract(i,"[0-9]{4}_[0-9]{2}_[0-9]{2}")
    ##
    s <- unlist(strsplit(sub(".fcs","",basename(i)),split = name.split))
    s <- grep("CONCATENATED",s,invert = T,value = T)
    if(!any(grepl("+",s,fixed = T))){s <- c(s,"PBMC")}
    ##
    vec.tmp <- stats::setNames(rep(NA,5),nm=c("subject","visit","condition","batch","cell.type"))
    if(length(s)==length(vec.tmp)){
      vec.tmp[] <- s
    }else if(length(s)==4){
      vec.tmp[!names(vec.tmp) %in% "batch"] <- s
      vec.tmp["batch"] <- batch
    }
    vec.tmp['batch.date'] <- batch.date
    vec.tmp[c('global.id','total.events')] <- flowCore::read.FCSheader(i, keyword = c("GLOBAL ID","$TOT"))[[1]]
    return(vec.tmp)
  })
  mdat <- data.frame(do.call(rbind,split.names.vec))
  mdat[,c('condition','batch')] <- lapply(mdat[,c('condition','batch')],factor)
  mdat$visit = factor(mdat$visit,levels = c("V4", "V6", "V7", "Adult"))
  mdat$total.events <- as.numeric(mdat$total.events)
  mdat$sample <- paste(mdat$subject, mdat$visit, mdat$condition, mdat$batch, mdat$cell.type, sep = "_")
  mdat$sample[grep("HD",mdat$sample,invert = T)] <- sub("_[0-9]{3}_","_",mdat$sample[grep("HD",mdat$sample,invert = T)])
  return(mdat)
}

#' List .fcs file paths
#'
#' @param root.dir.type c("modified","source")
#' @param sub.dir.type workflow directory c("normalized","trimmed","coded","debarcoded","concatenated","subsets")
#' @param subset.type subset population c("CD4+","CD8+","CD14+","CD19+","CD56+","TCRgd+")
#'
#' @return a batch-named list of .fcs file paths
#' @export
#'
get.fcs.paths.echo <- function(root.dir.type = c("modified","source"),sub.dir.type = NULL, subset.type = NULL){
  root <- match.arg(root.dir.type)
  fcs.dirs <- list.dirs(file.path(paste("data",root,sep = "_")))
  fcs.dirs <- grep(paste0(sub.dir.type,"$"),fcs.dirs,value = T)
  fcs.files <- sapply(fcs.dirs,function(i){
    fcs.files <- list.files(i,recursive=T,pattern=".fcs",full.names = T)
    fcs.files <- grep(subset.type,fcs.files,value = T)
  })
  fcs.files <- fcs.files[sapply(fcs.files,length)!=0]
  names(fcs.files) <- stringr::str_extract(names(fcs.files),"[0-9]{4}_[0-9]{2}_[0-9]{2}_ECHO_[0-9]{3}")
  return(fcs.files)
}

#' Get .fcs marker names from a list of .fcs files
#'
#' @param fcs.list list of .fcs files; returned from 'read.fcs.selected.markers'
#'
#' @return named character: metal_marker
#' @export
#'
get.metal.markers <- function(fcs.list){
  pn.names <- unique(lapply(fcs.list,flowCore::colnames));if(length(pn.names)==1) pn.names <- pn.names[[1]]
  ps.names <- unique(lapply(fcs.list,flowCore::markernames));if(length(ps.names)==1) ps.names <- ps.names[[1]]
  ##
  if(class(ps.names)=="list"){
    isometal <- unique(lapply(ps.names,function(i){
      tmp.string <- sapply(strsplit(i,"_"),'[[',1)
      tmp.string <- paste0(stringr::str_extract(tmp.string,"[0-9]{3}"),stringr::str_extract(tmp.string,"[A-Z][a-z]"))
    })
    )
    if(length(isometal)==1){
      isometal <- unlist(isometal)
      ps.names <- ps.names[unlist(lapply(ps.names,function(i) all(grepl(paste0(isometal,collapse = "|"),i))))][[1]]
    }
  }
  if(class(ps.names)=="character"&all(pn.names %in% names(ps.names))){
    return(ps.names)
  }
}

#' Generate a trim value for 'Time'
#'
#' @param time.vec numeric time vector
#' @param plot plot/visualize trim values
#'
#' @return a numeric value
#' @export
#'
#' @examples
time.trim <- function(time.vec,plot=F){
  #time histogram
  h <- hist(time.vec,breaks=200,plot=F)
  #low count bin cut; usually at end of acquisition
  cut.bins <- mean(h$counts)*.2
  if(all(!h$counts<cut.bins)){
    low.bin.break <-max(h$breaks)
  }else{
    low.bin.break <- which(h$counts<cut.bins)[1]
    #time at which low counts start
    low.bin.break <- h$breaks[low.bin.break]
  }
  ##
  trimmed.counts <- h$counts[h$breaks<low.bin.break]
  trimmed.counts.mean <- mean(trimmed.counts)
  trimmed.counts.break <- which(rev(trimmed.counts)<trimmed.counts.mean)[1]
  ##
  time.break <- rev(h$breaks[1:length(trimmed.counts)])[8]
  ##
  if(plot){
    plot(h)
    abline(h=trimmed.counts.mean,col="blue")
    abline(v=time.break,col="red")
  }
  else{
    return(time.break)
  }
}
