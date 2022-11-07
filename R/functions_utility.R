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

get.batches.fcs <- function(data.type.path, study.name){
  batches <- grep(study.name, list.dirs(data.type.path, recursive = F), value = T)
  batches <- sapply(batches, function(i) list.dirs(i, recursive = F), simplify = F)
  batches <- lapply(batches, function(i) list.files(i, pattern = ".fcs", full.names = TRUE, ignore.case = T))
  return(batches)
}

normalized.fcs_source.to.modified <- function(fcs.file.path, remove.file = FALSE){

  if(flowCore::read.FCSheader(fcs.file.path, keyword = "$CYTSN") == "StandAlone_Software"){

    print(paste("Moving/renaming:", fcs.file.path))

    dir.mod <- sub("source", "modified", dirname(fcs.file.path))
    dir.mod <- file.path(dirname(dir.mod), "01_normalized", basename(dir.mod))

    if(!dir.exists(dir.mod)){
      dir.create(dir.mod, recursive = T)
    }

    input.cytoffcs <- flowCore::read.FCS(fcs.file.path, transformation = F, truncate_max_range = F)
    input.cytoffcs@description$NORMALIZED <- "YES"
    input.cytoffcs@description$NORM.METHOD <- "EQ-P13H2302_ver2"

    FIL <- input.cytoffcs@description$`$FIL`

    batch.date.string.file <- stringr::str_extract(FIL,"[0-9]{8}")
    batch.date.string.dir <- gsub("_","",stringr::str_extract(dir.mod,"[0-9]{4}_[0-9]{2}_[0-9]{2}"))
    if(batch.date.string.file!=batch.date.string.dir){
      FIL <- sub("[0-9]{8}",batch.date.string.dir,FIL)
    }#needed to fix naming convention error;wrong date string

    FIL <- sub(";","",FIL)#needed to fix naming convention error
    FIL <- sub("E.*O", "ECHO", FIL)#needed to fix naming convention error
    #FIL <- sub("[0-9].fcs", "normalized.fcs", FIL)#drops CyTOF software appended number(s);new suffix
    #FIL <- gsub("_[0-9]{2}_", "_", FIL)#drops CyTOF software appended number(s)
    FIL <- sub("1_2", "2", FIL)#needed to fix naming convention error
    FIL <- sub("1_3_", "3", FIL)#needed to fix naming convention error
    FIL <- gsub(" - ", "_", FIL)#needed to fix naming convention error
    FIL <- gsub("-", "_", FIL)#needed to fix naming convention error
    FIL <- sub("_UNSTM_", "_UNSTIM_", FIL)#needed to fix naming convention error
    FIL <- sub("_Unstim_", "_UNSTIM_", FIL)#needed to fix naming convention error
    FIL <- sub("_STIM_", "_SEB_", FIL)#needed to fix naming convention error

    FIL <- stringr::str_extract(FIL, "ECHO_[0-9]{8}_[A-Z]+_[A-Z]{1}[0-9]{1}")
    FIL <- paste(FIL, "normalized.fcs", sep = "_")

    #aliquot.string <- stringr::str_extract(FIL, "_[A-Z]{1}[0-9]{1}_")
    aliquot.string <- substr(sub("ECHO_[0-9]{8}_[A-Z]+_", "", FIL),1,2)
    aliquot.letter <- stringr::str_extract(aliquot.string, "[A-Z]{1}")

    # if(any(grepl(paste0(aliquot.letter, "[0-9]{1}"), list.files(dir.mod, pattern = ".fcs")))){
    #   existing.aliquots <- grep(paste0(aliquot.letter, "[0-9]{1}"), list.files(dir.mod, pattern = ".fcs"), value = T)
    #   aliquot.strings <- sapply(existing.aliquots, FUN = stringr::str_extract, "_[A-Z]{1}[0-9]{1}_", USE.NAMES = F)
    #   aliquot.increment <- max(as.numeric(sapply(aliquot.strings, FUN = stringr::str_extract, "[0-9]{1}", USE.NAMES = F)))+1
    #   aliquot.string <- sub("[0-9]{1}", aliquot.increment, aliquot.string)
    #   FIL <- sub("_[A-Z]{1}[0-9]{1}_", aliquot.string, FIL)
    # }

    if(any(grepl(aliquot.string, list.files(dir.mod, pattern = ".fcs")))){
      existing.aliquots <- grep(paste0(aliquot.letter, "[0-9]{1}"), list.files(dir.mod, pattern = ".fcs"), value = T)
      aliquot.strings <- sapply(existing.aliquots, FUN = stringr::str_extract, "_[A-Z]{1}[0-9]{1}_", USE.NAMES = F)
      aliquot.increment <- max(as.numeric(sapply(aliquot.strings, FUN = stringr::str_extract, "[0-9]{1}", USE.NAMES = F)))+1
      aliquot.string <- sub("[0-9]{1}", aliquot.increment, aliquot.string)
      FIL <- sub("[A-Z]{1}[0-9]{1}", aliquot.string, FIL)
    }

    input.cytoffcs@description$`$FIL` <- FIL

    if(grepl("-", input.cytoffcs@description$GUID)){
      GUID <- sub("-", "_", gsub(" ", "", input.cytoffcs@description$GUID))
    }else{
      GUID <- input.cytoffcs@description$GUID
    }

    FIL <- input.cytoffcs@description$`$FIL`

    STUDY <- unlist(strsplit(FIL, "_"))[1]#
    POOL <- unlist(strsplit(FIL, "_"))[2]#
    CONDITION <- unlist(strsplit(FIL, "_"))[3]#
    if(CONDITION == "STIM"){
      CONDITION <- "SEB"#needed to fix naming convention error
    }
    CONDITION <- toupper(CONDITION)#needed to fix naming convention error

    if(length(unlist(strsplit(FIL, "_"))) >= 6){
      POOL.ALIQUOT <- paste(unlist(strsplit(FIL, "_"))[4:5], collapse = "")#needed to fix naming convention error
    }else if(length(unlist(strsplit(FIL, "_"))) == 5){
      POOL.ALIQUOT <- unlist(strsplit(FIL, "_"))[4]#
    }

    ALIQUOT.TOTAL <- nrow(input.cytoffcs)

    input.cytoffcs@description$STUDY <- STUDY
    input.cytoffcs@description$POOL <- POOL
    input.cytoffcs@description$CONDITION <- CONDITION
    input.cytoffcs@description$POOL.ALIQUOT <- POOL.ALIQUOT
    input.cytoffcs@description$ALIQUOT.TOTAL <- ALIQUOT.TOTAL
    input.cytoffcs@description$EXP <- sub("ECHO_", "", stringr::str_extract(fcs.file.path, "ECHO_[0-9]{3}"))

    flowCore::write.FCS(input.cytoffcs, filename = file.path(dir.mod, input.cytoffcs@description$`$FIL`))

    new.file.exists <- file.exists(file.path(dir.mod, input.cytoffcs@description$`$FIL`))
    new.file.same.total.as.old.file <- unlist(flowCore::read.FCSheader(fcs.file.path, keyword = "$TOT")) ==
      unlist(flowCore::read.FCSheader(file.path(dir.mod, input.cytoffcs@description$`$FIL`), keyword = "$TOT"))

    if(new.file.exists & new.file.same.total.as.old.file & remove.file == TRUE){
      message(paste("Removing:", fcs.file.path))
      file.remove(fcs.file.path)
    }
  }else{
    print("Source data...")
  }

  # sapply(unlist(echo.fcs.to.move), function(i) {
  #
  #   dir.mod <- sub("source", "modified", dirname(i))
  #   dir.mod <- file.path(dirname(dir.mod), "01_normalized", basename(dir.mod))
  #
  #   if(!dir.exists(dir.mod)){
  #     dir.create(dir.mod, recursive = T)
  #   }
  #
  #   input.cytoffcs <- read.FCS(i, transformation = F, truncate_max_range = F)
  #   input.cytoffcs@description$NORMALIZED <- "YES"
  #   input.cytoffcs@description$NORM.METHOD <- "EQ-P13H2302_ver2"
  #   input.cytoffcs@description$`$FIL` <- sub(paste0(paste0("_", c(0,1,2)), ".fcs", collapse = "|"),
  #                                            "_normalized.fcs", input.cytoffcs@description$`$FIL`)#drops CyTOF software appended number
  #   input.cytoffcs@description$`$FIL` <- sub("1_2", "2", input.cytoffcs@description$`$FIL`)#needed to fix naming convention error
  #   input.cytoffcs@description$`$FIL` <- sub("1_3", "3", input.cytoffcs@description$`$FIL`)#needed to fix naming convention error
  #   input.cytoffcs@description$`$FIL` <- sub(" - ", "_", input.cytoffcs@description$`$FIL`)#needed to fix naming convention error
  #   input.cytoffcs@description$`$FIL` <- sub("_Unstim_", "_UNSTIM_", input.cytoffcs@description$`$FIL`)#needed to fix naming convention error
  #   input.cytoffcs@description$`$FIL` <- sub("_STIM_", "_SEB_", input.cytoffcs@description$`$FIL`)#needed to fix naming convention error
  #
  #   if(grepl("-", input.cytoffcs@description$GUID)){
  #     GUID <- sub("-", "_", gsub(" ", "", input.cytoffcs@description$GUID))
  #   }else{
  #     GUID <- input.cytoffcs@description$GUID
  #   }
  #
  #   FIL <- input.cytoffcs@description$`$FIL`
  #
  #   STUDY <- unlist(strsplit(FIL, "_"))[1]#
  #   POOL <- unlist(strsplit(FIL, "_"))[2]#
  #   CONDITION <- unlist(strsplit(FIL, "_"))[3]#
  #   if(CONDITION == "STIM"){
  #     CONDITION <- "SEB"#needed to fix naming convention error
  #   }
  #   CONDITION <- toupper(CONDITION)#needed to fix naming convention error
  #
  #   if(length(unlist(strsplit(FIL, "_"))) >= 6){
  #     POOL.ALIQUOT <- paste(unlist(strsplit(FIL, "_"))[4:5], collapse = "")#needed to fix naming convention error
  #   }else if(length(unlist(strsplit(FIL, "_"))) == 5){
  #     POOL.ALIQUOT <- unlist(strsplit(FIL, "_"))[4]#
  #   }
  #
  #   ALIQUOT.TOTAL <- nrow(input.cytoffcs)
  #
  #   input.cytoffcs@description$STUDY <- STUDY
  #   input.cytoffcs@description$POOL <- POOL
  #   input.cytoffcs@description$CONDITION <- CONDITION
  #   input.cytoffcs@description$POOL.ALIQUOT <- POOL.ALIQUOT
  #   input.cytoffcs@description$ALIQUOT.TOTAL <- ALIQUOT.TOTAL
  #   input.cytoffcs@description$EXP <- exp.number
  #
  #   write.FCS(input.cytoffcs, filename = file.path(dir.mod, input.cytoffcs@description$`$FIL`))
  #
  #   #file.remove(i)
  # })
}

fcs.headers.no.pars <- function(fcs.file.paths, keywords.build = TRUE, return.combined = TRUE){
  fcs.header <- sapply(fcs.file.paths, function(i){
    i <- as.data.frame(flowCore::read.FCSheader(i))
    i <- data.frame(keyword = rownames(i),
                    value = i[,1])
    i <- i[-grep("\\$P[0-9]{1}", i$keyword), ]
    i <- i[-grep("\\$BEGIN", i$keyword), ]
    i <- i[-grep("\\$BYTEORD", i$keyword), ]
    i <- i[-grep("\\$DATATYPE", i$keyword), ]
    i <- i[-grep("\\$END", i$keyword), ]
    i <- i[-grep("\\$MODE", i$keyword), ]
    i <- i[-grep("\\$NEXTDATA", i$keyword), ]
    i <- i[-grep("\\$PAR", i$keyword), ]

    helios.keywords <- c("FluidigmBarcoded",
                         "FluidigmMassTagsChecked",
                         "FluidigmNormalized",
                         "TargetGemStoneMethod")

    if(any(helios.keywords %in% i$keyword)){
      i <- i[-grep(paste0(helios.keywords, collapse = "|"), i$keyword), ]
    }

    i
  }, USE.NAMES = T, simplify = F)

  if(unique(sapply(fcs.header, nrow)) == length(Reduce(union, sapply(fcs.header, function(i) i$keyword)))&keywords.build == TRUE){
    keywords.build <- c("$BTIM","$CYT","$CYTSN","$DATE","$ETIM","$FIL","$TOT","ALIQUOT.TOTAL","CONDITION","EXP","FCSversion",
                        "FILENAME","GUID","NORM.METHOD","NORMALIZED","ORIGINALGUID","POOL","POOL.ALIQUOT","STUDY")

    fcs.header <- sapply(fcs.header, function(i){
      if(!all(keywords.build %in% i$keyword)){
        stop("Keyword conflict...")
      }else{
        i <- i[i$keyword %in% keywords.build, ]
      }
    }, simplify = F)
  }

  if(return.combined == TRUE){
    keywords.combined <- stats::setNames(vector(mode = "character", length = unique(sapply(fcs.header, nrow))),
                                         nm = unlist(unique(lapply(fcs.header, function(i) i$keyword))))
    for(k in seq(keywords.combined)){
      keywords.combined[k] <- paste0(sapply(fcs.header, function(i) i$value[i$keyword == names(keywords.combined)[k]]), collapse = ";")
    }

    keywords.combined <- sapply(keywords.combined, function(i){
      if(length(unique(strsplit(i, ";")[[1]]))==1){
        i <- unique(strsplit(i, ";")[[1]])
      }else{
        i <- i
      }
    })
    pluralize.names <- names(keywords.combined)[sapply(keywords.combined, function(i) grepl(";", i))]
    names(keywords.combined)[names(keywords.combined) %in% pluralize.names] <- paste0(sub("\\$", "", pluralize.names), "S")
    keywords.combined
  }else{
    fcs.header
  }
}

iso.metal.marker.fix <- function(marker.vec,name.split="_"){
  m.split <- strsplit(marker.vec,name.split);m.split1 <- sapply(m.split,'[[',1);m.split2 <- sapply(m.split,'[[',2)
  metal <- stringr::str_extract(m.split1,"[A-Za-z]+")
  iso <- stringr::str_extract(m.split1,"[0-9]+")
  iso.metal.marker <- paste(paste0(iso,metal),m.split2,sep = "_")
  iso.metal.marker <- factor(iso.metal.marker)
  return(iso.metal.marker)
}

get.range.limits <- function(fcs.filepaths,marker.only.names=F,cofactor=5,markers.vec=NULL){
  ranges<-sapply(fcs.filepaths,function(f.path){
    header <- flowCore::read.FCSheader(f.path)[[1]]
    stains <- grep("_",header[grep("P[0-9]+S",names(header),value = T)],value = T)
    ranges <- stats::setNames(as.numeric(header[sub("S","R",names(stains))]),
                              nm=stains
    )
    names(ranges) <- iso.metal.marker.fix(names(ranges))
    if(!any(grepl("151Eu_CD69",names(ranges)))){
      ranges <- c(ranges,"151Eu_CD69"=0)
    }
    if(!any(grepl("155Gd_TCRgd",names(ranges)))){
      ranges <- c(ranges,"155Gd_TCRgd"=0)
    }
    ranges<-ranges[sort(names(ranges))]
    marker.names <- sapply(strsplit(names(ranges),"_"),'[[',2)
    if(marker.only.names){
      names(ranges) <- marker.names
    }
    if(!is.null(markers.vec)){
      ranges<-ranges[marker.names%in%markers.vec]
    }
    return(ranges)
  },simplify = F)
  ranges <- do.call(rbind,ranges)
  range.limits <- apply(ranges,2,function(x,cf=cofactor){asinh(mean(x/cf))})
  return(range.limits)
}
