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

#' Get .fcs marker names from .fcs file path(s); attempts to be agnostic to naming convention
#'
#' @param fcs.path .fcs file path
#' @param selected.markers optional character vector of selected marker names; for getting '$P##N' column names
#' @param suppress.check logical; set to TRUE to suppress error checking and return incomplete matching of selected.markers
#'
#' @return a list; marker 'names' and marker 'stains';stops if marker not found
#' @export
#'
fcs.markers.agnostic <- function(fcs.path,selected.markers=NULL,suppress.check=F){
  header <- flowCore::read.FCSheader(fcs.path)[[1]]
  p <- sapply(c("N","S"),function(i){
    p <- header[grep(paste0("P[0-9]+",i),names(header),value = T)]
    p <- p[order(as.numeric(stringr::str_extract(names(p),"[0-9]+")))]
  })

  pn.selected <- p$N[which(p$N %in% selected.markers)]

  common.split.counts <- sapply(c(".","_","-"," "),function(split) sum(grepl(split,p$S,fixed = T)))
  if(max(common.split.counts)!=length(p$S)){
    most.likely.split <- ""
  }else{
    most.likely.split <- names(common.split.counts)[which.max(common.split.counts)]
  }

  if(!is.null(selected.markers)){
    ps.selected <- names(grep(paste0(most.likely.split,paste0(selected.markers,"$"),collapse = "|"),p$S,value = T))
    if(length(ps.selected)==0){
      ps.selected <- names(grep(paste(paste0(selected.markers,most.likely.split),collapse = "|"),p$S,value = T))
      if(length(ps.selected)==0){
        stop("Selected marker(s) not found")
      }
    }
    pn.selected <- c(p$N[names(p$N) %in% sub("S","N",ps.selected)],pn.selected)
    if(length(selected.markers)!=length(pn.selected)&!suppress.check){
      stop(paste("Selected",length(selected.markers), "markers but found:",length(pn.selected)))
    }else{
      p$N <- p$N[p$N %in% pn.selected]
      p$S <- p$S[names(p$S) %in% ps.selected]
      return(p)
    }
  }else{
    return(p)
  }
}

#' Check if markers exist in all .fcs files
#'
#' @param fcs.paths .fcs file paths
#' @param markers.vec a character vector of marker names
#' @param ... additional arguments to pass to fcs.markers
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

#' Check if selected markers exist in all .fcs files
#'
#' @param fcs.paths .fcs file paths
#' @param selected.markers a character vector of select marker names
#'
#' @return logical test; stops if a non-unique list of markers/missing markers are found
#'
fcs.markers.agnostic.check <- function(fcs.paths,selected.markers){
  markers.list <- lapply(fcs.paths, function(i) {
    fcs.markers.agnostic(i,selected.markers)
  })
  unique.n <- unique(lapply(markers.list,'[[','N'))
  unique.n.intersect <- Reduce(intersect,unique.n)
  if(!is.vector(unique.n.intersect)){
    stop("Unique names ('P$N') error: need a vector...")
  }
  ##
  unique.s <- unique(lapply(markers.list,'[[','S'))
  unique.s.length <- unique(sapply(unique.s,length))
  ##
  if(length(unique.n.intersect)!=unique.s.length){
    print(unique.n)
    print(unique.s)
    print(unique.n.intersect)
    print(unique.s.length)
    stop(paste("Length of unique names ('P$N') does not equal length of stains ('P$S'))",
               "The provided selected markers are not common/found across this list of .fcs files",
               sep = "\n")
    )
  }
}

#' Read .fcs file; keep only selected markers
#'
#' @param fcs.path .fcs file path
#' @param selected.markers character vector of marker names to keep
#' @param include.scatter logical; for fluor/flow based assays, include light-scatter channels
#' @param trim.time logical; for fluor/flow based assays, trim time from end
#' @param trim.scatter logical; for fluor/flow based assays, trim high light-scatter values
#' @param comp logical; for fluor/flow based assays, compensate; sample-specific
#' @param comp.mat logical; for fluor/flow based assays, compensate using a provided (modified) matrix
#' @param subsample.val numeric value; if provided, will subsample the fcs file
#' @param suppress.check logical; set to TRUE to suppress error checking and return incomplete matching of selected.markers
#'
#' @return a flowFrame; no transformation applied; no truncation; pre-processed according to arguments
#' @export
#'
read.fcs.selected.markers <- function(fcs.path, selected.markers,include.scatter=T,trim.time=T,trim.scatter=T,comp=T,comp.mat=NULL,subsample.val=NULL,suppress.check=F){
  selected.markers.names <- fcs.markers.agnostic(fcs.path,selected.markers=selected.markers,suppress.check = suppress.check)$N#return the 'names' from the list
  if(grepl("LASER",flowCore::read.FCSheader(fcs.path))){
    if(include.scatter){
      selected.markers.names <- c(selected.markers.names,
                                  paste(rep(c('FSC','SSC'),each=3),c('A','H','W'),sep = "-")
      )
    }
    fcs <- flowCore::read.FCS(fcs.path,transformation=F,truncate_max_range=F)
    if(trim.time){
      fcs@exprs <- fcs@exprs[fcs@exprs[,'Time']<trim.time.value(fcs@exprs[,'Time']),]
    }
    if(trim.scatter){
      scatter.trim.list <- sapply(paste(rep(c('FSC','SSC'),each=3),c('A','H','W'),sep = "-"),function(scatter){
        which(fcs@exprs[,scatter]>250000)
      })
      fcs@exprs <- fcs@exprs[-Reduce(union,scatter.trim.list),]
    }
    if(comp&is.null(comp.mat)){
      if(!is.null(fcs@description$SPILL)){
        fcs <- flowCore::compensate(fcs,fcs@description$SPILL)
      }else if(!is.null(fcs@description$`$SPILLOVER`)){
        fcs <- flowCore::compensate(fcs,fcs@description$`$SPILLOVER`)
      }
    }else if(comp&!is.null(comp.mat)){
      fcs <- flowCore::compensate(fcs,comp.mat)
    }
    fcs <- fcs[,flowCore::colnames(fcs) %in% selected.markers.names]
  }else{
    fcs <- flowCore::read.FCS(fcs.path,
                              column.pattern = paste0(selected.markers.names,collapse = "|"),
                              transformation = F, truncate_max_range = F)
  }
  ##
  if(!is.null(subsample.val)){
    if(nrow(fcs)>subsample.val){
      set.seed(20040501)
      row.samp <- sample(1:nrow(fcs), subsample.val)
      fcs <- fcs[row.samp,]
    }
  }
  ##
  return(fcs)
}

#' Read .fcs files; keep only selected markers; parallelized
#'
#' @param fcs.paths .fcs file paths
#' @param ... additional arguments as handled by read.fcs.selected.markers
#' @param list.check logical; set to FALSE to suppress 'fcs.markers.agnostic.check'
#'
#' @return a list of individual flowFrames
#' @export
#'
read.fcs.selected.markers.parallel <- function(fcs.paths,list.check=T,...){
  #make a check before initializing parallel clusters
  args <- list(...)
  if("selected.markers" %in% names(args)){
    selected.markers <- args$selected.markers
  }
  if(list.check){
    fcs.markers.agnostic.check(fcs.paths,selected.markers)
  }
  #
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, envir = environment(), c("..."))
  fcs.list <- parallel::parSapply(cl, fcs.paths, function(i,...){
    fcs <- read.fcs.selected.markers(fcs.path=i,...)
  },...)
  parallel::stopCluster(cl)
  return(fcs.list)
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
  if(!is.character(fcs.paths)&!is.vector(fcs.paths)){
    stop("fcs.paths needs to be a character vector of .fcs file paths")
  }
  split.names.vec <- lapply(fcs.paths, function(i){
    batch <- sub("ECHO_","",stringr::str_extract(i,"ECHO_[0-9]{3}"))
    batch.date <- stringr::str_extract(i,"[0-9]{4}_[0-9]{2}_[0-9]{2}")
    ##
    s <- unlist(strsplit(sub(".fcs","",basename(i)),split = name.split))
    s <- grep("CONCATENATED",s,invert = T,value = T)
    ##
    if(length(s)>grep("\\+|PBMC",s)){
      s <- s[1:grep("\\+|PBMC",s)]
    }
    ##
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
    if(!is.null(subset.type)){
      fcs.files <- grep(subset.type,fcs.files,value = T)
    }
    return(fcs.files)
  },simplify = F)
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
  if(is.list(ps.names)){
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
  if(is.character(ps.names)&all(pn.names %in% names(ps.names))){
    return(ps.names)
  }
}

#' Generate a trim value for 'Time'
#'
#' @param time.vec numeric time vector
#' @param plot plot/visualize trim value; histogram
#'
#' @return a numeric value
#'
trim.time.value <- function(time.vec,plot=F){
  #time histogram
  h <- graphics::hist(time.vec,breaks=200,plot=F)
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
    graphics::abline(h=trimmed.counts.mean,col="blue")
    graphics::abline(v=time.break,col="red")
  }
  else{
    return(time.break)
  }
}
