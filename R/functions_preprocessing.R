bead.viability.trim <- function(normalized.fcs.file.path, viability_metal, bead_metal, cut.method = "counts"){

  dir.mod <- sub("01_normalized", "02_trimmed", dirname(normalized.fcs.file.path))

  if(!dir.exists(dir.mod)){
    dir.create(dir.mod, recursive = T)
  }

  dir.mod.p <- file.path(dirname(dir.mod), "trim_plots")

  if(!dir.exists(dir.mod.p)){
    dir.create(dir.mod.p, recursive = T)
  }

  input.cytoffcs <- flowCore::read.FCS(normalized.fcs.file.path, transformation = F, truncate_max_range = F)

  viability.beads <- stats::setNames(c(paste0(viability_metal, "Di"), paste0(bead_metal, "Di")), c("viability", "beads"))

  fcs.tmp <- input.cytoffcs[, which(flowCore::colnames(input.cytoffcs) %in% c(viability.beads, "Time"))]
  fcs.tmp <- flowCore::transform(fcs.tmp, flowCore::transformList(viability.beads, list(function(x, cofactor = 5){ asinh(x)/cofactor})))

  #
  if(cut.method == "densities"){
    bv.densities <- list(b = stats::density(subset(fcs.tmp@exprs[, viability.beads["beads"]],
                                            fcs.tmp@exprs[, viability.beads["beads"]] > stats::quantile(fcs.tmp@exprs[, viability.beads["beads"]], probs = .95))),
                         v = stats::density(subset(fcs.tmp@exprs[, viability.beads["viability"]],
                                            fcs.tmp@exprs[, viability.beads["viability"]] > 0))
    )

    valley <- function(v){
      #v.sub <- v$y[which.max(v$y):length(v$y)]#find first highest peak and subset
      v.sub <- v$y[(which((v$y == cummax(v$y)) == FALSE)[1]-1):length(v$y)]#find first peak and subset
      print(max(v.sub))

      valley.y.min <- v.sub[which((v.sub == cummin(v.sub)) == FALSE)[1]-1]#tests for monotonic decrease, finding trough/valley value

      valley.x <- v$x[which(v$y %in% valley.y.min)]#find value at trough/valley position
      valley.x
    }

    trim.values <- lapply(bv.densities, valley)

    viability.cutoff <- stats::quantile(fcs.tmp@exprs[, viability.beads["viability"]][which(fcs.tmp@exprs[, viability.beads["viability"]] < trim.values$v)],
                                 probs = c(0.98))

    bead.cutoff <- stats::quantile(fcs.tmp@exprs[, viability.beads["beads"]][which(fcs.tmp@exprs[, viability.beads["beads"]] < trim.values$b)],
                            probs = c(0.99))

    trim.index <- intersect(which(fcs.tmp@exprs[, viability.beads["beads"]] < bead.cutoff),
                            which(fcs.tmp@exprs[, viability.beads["viability"]] < viability.cutoff)
    )
    trim.index
  }else if(cut.method == "counts"){
    peak.detection.histogram.counts <- function(dat.vector, cut.percent = 0.65){
      x <- graphics::hist(dat.vector*-1, breaks = 200, plot = FALSE)
      peak.cut <- c(25,50,100,200,500,1000,2000,5000)
      peak <- vector(mode = "numeric", length = length(peak.cut))
      for(i in seq(peak)){
        peak[[i]] <- x$mids[x$counts>peak.cut[i]][which((x$counts[x$counts>peak.cut[i]] == cummax(x$counts[x$counts>peak.cut[i]])) == FALSE)[1]-1]*-1
      }
      peak <- peak[peak>1]
      print(peak)
      peak <- as.numeric(names(sort(table(peak),decreasing=TRUE)[1]))*cut.percent
      peak
    }

    bead.cutoff <- peak.detection.histogram.counts(fcs.tmp@exprs[fcs.tmp@exprs[, viability.beads["beads"]]>0, viability.beads["beads"]])

    viability.cutoff <- peak.detection.histogram.counts(fcs.tmp@exprs[fcs.tmp@exprs[, viability.beads["viability"]]>0, viability.beads["viability"]])

    trim.index <- intersect(which(fcs.tmp@exprs[, viability.beads["beads"]] < bead.cutoff),
                            which(fcs.tmp@exprs[, viability.beads["viability"]] < viability.cutoff)
    )
    trim.index
  }else{
    stop("'cut.method must be either 'counts' or 'densities'")
  }

  input.cytoffcs@exprs <- input.cytoffcs@exprs[trim.index, ]
  ##
  input.cytoffcs@description$`$FIL` <- sub("normalized.fcs", "normalized_trimmed_mod2.fcs", input.cytoffcs@description$`$FIL`)

  flowCore::write.FCS(input.cytoffcs, file.path(dir.mod, input.cytoffcs@description$`$FIL`))
  ##
  dats <- list(dat.og = data.frame(fcs.tmp@exprs, type = "original"),
               dat.mod = data.frame(fcs.tmp@exprs[trim.index, ], type = "trimmed")
  )

  dats <- lapply(dats, function(i) {
    colnames(i)[which(colnames(i) %in% viability.beads)] <- names(viability.beads[match(colnames(i)[which(colnames(i) %in% viability.beads)],
                                                                                        viability.beads)])
    i
  })

  event.total <- data.frame(total = sapply(dats, nrow), type = c("original", "trimmed"))
  dats <- data.table::rbindlist(dats)

  plot.viability <- ggplot2::ggplot(dats[dats$viability>0,], ggplot2::aes(!!quote(Time), !!quote(viability))) +
    ggplot2::geom_hex() +
    viridis::scale_fill_viridis(option = "plasma", limits = c(100, 4000), oob = scales::squish) +
    ggplot2::geom_hline(yintercept = viability.cutoff, color = "red") +
    #geom_hline(yintercept = viability.peak.value, color = "red") +
    ggplot2::facet_wrap(~type) +
    ggplot2::geom_text(data = event.total, ggplot2::aes(x = Inf, y = Inf, label = !!quote(total)), hjust = 1.5, vjust = 1.5) +
    ggplot2::labs(title = basename(normalized.fcs.file.path))

  plot.beads <- ggplot2::ggplot(dats[dats$beads>0,], ggplot2::aes(!!quote(Time), !!quote(beads))) +
    ggplot2::geom_hex() +
    viridis::scale_fill_viridis(option = "plasma", limits = c(50, 500), oob = scales::squish) +
    ggplot2::geom_hline(yintercept = bead.cutoff, color = "red") +
    #geom_hline(yintercept = bead.peak.value, color = "red") +
    ggplot2::facet_wrap(~type) +
    ggplot2::geom_text(data = event.total, ggplot2::aes(x = Inf, y = Inf, label = !!quote(total)), hjust = 1.5, vjust = 1.5) +
    ggplot2::labs(title = basename(normalized.fcs.file.path))

  plot.grobs <- gridExtra::arrangeGrob(plot.viability, plot.beads)
  ggplot2::ggsave(file = file.path(dir.mod.p, sub("_normalized.fcs", ".pdf", basename(normalized.fcs.file.path))),
                  plot.grobs, width = 9, height = 7)
}

combine.trim.plots <- function(echo.exp.number){
  exp.dir <- grep(echo.exp.number, list.dirs("../Mass_Cytometry/data_modified/", recursive = F), value = T)
  trim.plot.dir <- grep("trim_plots", list.dirs(exp.dir, recursive = T), value = T)

  out.file <- file.path(trim.plot.dir, paste(basename(exp.dir), "combined_trim.pdf", sep = "_"))
  if(file.exists(out.file)){
    file.remove(out.file)
  }

  trim.plots <- list.files(trim.plot.dir, pattern = ".pdf", full.names = T)

  pdftools::pdf_combine(trim.plots, output = out.file)
}

concatenate.fcs <- function(debarcoded.fcs.file.paths){
  debarcoded.fcs.file.paths <- grep("MISMATCH", debarcoded.fcs.file.paths, invert = T, value = T)#drop 'MISMATCH' files
  assigned.names <- stats::setNames(nm = basename(debarcoded.fcs.file.paths))
  name.split <- sapply(sub(".fcs", "", assigned.names), function(i) strsplit(i, "_"))
  pid.visit.condition <- unique(sapply(names(name.split), function(i){
    if(grepl("UNASSIGNED|Adult", i)){
      #"UNASSIGNED"
      paste(name.split[[i]][1],
            name.split[[i]][2],
            name.split[[i]][3],
            name.split[[i]][4],
            sep = "_")
    }else{
      paste(name.split[[i]][1],
            name.split[[i]][2],
            name.split[[i]][3],
            sep = "_")
    }
  })
  )

  dir.mod <- sub("04_debarcoded", "05_concatenated", unique(dirname(debarcoded.fcs.file.paths)))

  invisible(
    sapply(dir.mod, function(i){
      if(!dir.exists(i)){
        dir.create(i, recursive = T)
      }
    })
  )

  invisible(sapply(pid.visit.condition, function(i){

    fcs.paths.tmp <- grep(i, debarcoded.fcs.file.paths, value = T)

    keyword.headers <- fcs.headers.no.pars(fcs.paths.tmp, keywords.build = FALSE)

    par.num <- sapply(fcs.paths.tmp, function(i) flowCore::read.FCSheader(i, keyword = "$PAR")[[1]])

    if(length(unique(par.num))>1 & max(as.numeric(par.num)) - min(as.numeric(par.num)) == 1){
      message("Found 1 single parameter difference between one frame and the others...")
      affected.frame <- names(which(sapply(fcs.paths.tmp, FUN = flowCore::read.FCSheader, keyword = "$PAR", USE.NAMES = F) == names(table(par.num))[min(table(par.num))]))
      unaffected.frames <- names(which(sapply(fcs.paths.tmp, FUN = flowCore::read.FCSheader, keyword = "$PAR", USE.NAMES = F) == max(names(table(par.num)))))
      message(paste("Parameter issue with:", affected.frame, sep = " "))

      affected.PS <- flowCore::read.FCSheader(affected.frame, keyword = grep("\\$P.*S", names(flowCore::read.FCSheader(affected.frame)[[1]]), value = T))
      unaffected.PS <-flowCore::read.FCSheader(unaffected.frames[1], keyword = grep("\\$P.*S", names(flowCore::read.FCSheader(unaffected.frames[1])[[1]]), value = T))

      extra.par.index <- which(!unlist(unaffected.PS, use.names = F) %in% unlist(affected.PS, use.names = F))
      extra.par <- unlist(unaffected.PS, use.names = F)[extra.par.index]
      message(paste("Extra parameter in the other frames:", extra.par, sep = " "))

      col.pattern <- as.character(flowCore::read.FCSheader(unaffected.frames[1], keyword = sub("S", "N", names(unaffected.PS[[1]])[extra.par.index]))[[1]])
      message(paste("Reading in as a flowset but dropping the following:", col.pattern, sep = " "))
      input.cytofset <- flowCore::read.flowSet(fcs.paths.tmp, transformation = F, truncate_max_range = F, column.pattern = col.pattern, invert.pattern = T)
    }else if(length(unique(par.num))>1){
      message("Parameter number mismatch; reading in as multiple flowsets")
      input.cytofset <- vector(mode = "list", length = length(unique(par.num)))
      for(j in seq(input.cytofset)){
        input.cytofset[[j]] <- flowCore::read.flowSet(fcs.paths.tmp[par.num == unique(par.num)[[j]]], transformation = F, truncate_max_range = F)
      }
      message(paste("Dropping the following columns:",paste0(Reduce(setdiff, sapply(input.cytofset, colnames)),collapse = ","),sep = " "))
      shared.colnames <- Reduce(intersect, sapply(input.cytofset, flowCore::colnames))
      input.cytofset <- sapply(input.cytofset, function(i) i <- i[,shared.colnames])

      if(length(input.cytofset) <= 2){
        input.cytofset <- do.call(flowCore::rbind2, input.cytofset)
      }
    }else{
      input.cytofset <- flowCore::read.flowSet(fcs.paths.tmp, transformation = F, truncate_max_range = F)
      #shared.colnames <- Reduce(intersect, fsApply(input.cytofset, colnames, simplify = F))
      #input.cytofset <- fsApply(input.cytofset, function(f) f <- f[,shared.colnames])
    }

    if(suppressWarnings(is.list(flowCore::markernames(input.cytofset)))){
      message("Name-fixing parameters")
      input.cytofset <- flowCore::fsApply(input.cytofset, function(f){
        beads.index <- grep("beads", f@parameters@data$desc, ignore.case = T)
        f@parameters@data$desc[beads.index] <- sub("_.*","_Norm_beads",f@parameters@data$desc[beads.index])
        #
        f@parameters@data$desc <- sub("viablitiy", "Viability", f@parameters@data$desc)#name fix
        f@parameters@data$desc <- sub("viability", "Viability", f@parameters@data$desc)#name fix
        viability.index <- grep("viability", f@parameters@data$desc, ignore.case = T)
        if(grepl("Pt",f@parameters@data$desc[viability.index])){
          f@parameters@data$desc[viability.index] <- sub("_.*","_viability_cisplatin",f@parameters@data$desc[viability.index])
        }#name fix
        #
        f@parameters@data$desc <- sub("CD49a", "CD49d", f@parameters@data$desc)#name fix
        f@parameters@data$desc <- sub("TGFb1", "TGFbeta1", f@parameters@data$desc)#name fix
        f@parameters@data$desc <- sub("background", "Background", f@parameters@data$desc)#name fix
        f@parameters@data$desc <- sub("176Yb$", "176Yb_CD127", f@parameters@data$desc)#name fix
        f@parameters@data$desc <- sub("209Bi$", "209Bi_CD57", f@parameters@data$desc)#name fix
        f@parameters@data$desc[grep("190BCKG",f@parameters@data$desc)] <- "190BCKG_Noise"#name fix

        desc.na.index <- !f@parameters@data$name %in% c("Time", "Event_length", "node_CD45", "barcode", "barcode_cut") & is.na(f@parameters@data$desc)
        f@parameters@data$desc[desc.na.index] <- f@parameters@data$name[desc.na.index]
        f@parameters@data$desc <- sub("Di", "", f@parameters@data$desc)#name fix

        marker.index <- grepl("_", f@parameters@data$desc)
        marker.metal <- sapply(strsplit(f@parameters@data$desc[marker.index], "_"), '[', 1)
        marker.metal <- paste0(gsub("[A-Za-z]", "",marker.metal),
                               gsub("[0-9]", "", marker.metal))
        marker.marker <- sapply(strsplit(f@parameters@data$desc[marker.index], "_"), '[', 2)
        f@parameters@data$desc[marker.index] <- paste(marker.metal, marker.marker, sep = "_")

        metal.index <- !grepl("_", f@parameters@data$desc) & !is.na(f@parameters@data$desc)
        f@parameters@data$desc[metal.index] <- paste0(gsub("[A-Za-z]", "", f@parameters@data$desc[metal.index]),
                                                      gsub("[0-9]", "", f@parameters@data$desc[metal.index]))

        f
      })
    }

    tmp.dat <- flowCore::fsApply(input.cytofset, function(f){
      dat <- f@parameters@data[,c("name", "desc")]
      rownames(dat) <- NULL
      dat
    })

    if(all(sapply(sapply(tmp.dat, function(i) unname(i$name), simplify = F), FUN = identical, unname(tmp.dat[[length(tmp.dat)]]$name))) &
       all(sapply(sapply(tmp.dat, function(i) unname(i$desc), simplify = F), FUN = identical, unname(tmp.dat[[length(tmp.dat)]]$desc)))){
      tmp.dat <- tmp.dat[[1]]
    }else if(!all(sapply(tmp.dat, FUN = identical, tmp.dat[[length(tmp.dat)]]))){
      message(paste("Name or desc issue with:", names(which(!sapply(tmp.dat, FUN = identical, tmp.dat[[length(tmp.dat)]])))))
      affected.markernames <- flowCore::markernames(input.cytofset[[which(!sapply(tmp.dat, FUN = identical, tmp.dat[[length(tmp.dat)]]))]])
      unaffected.markernames <- flowCore::markernames(input.cytofset[which(sapply(tmp.dat, FUN = identical, tmp.dat[[length(tmp.dat)]]))])
      fix.these <- which(affected.markernames != unaffected.markernames)
      message(paste("These need fixing:", paste0(names(fix.these), collapse = ","), sep = " "))
      message(paste("The other frames have markernames:", paste0(unaffected.markernames[fix.these], collapse = ","), sep = " "))
      message("Dropping affected markernames; using correctly named markers in the concatenate...")
      tmp.dat <- tmp.dat[-which(!sapply(tmp.dat, FUN = identical, tmp.dat[[length(tmp.dat)]]))]
      message("Retesting for unified names/markers")
      if(!all(sapply(tmp.dat, FUN = identical, tmp.dat[[length(tmp.dat)]]))){
        stop("Either name or desc is not unified...")
      }else{
        tmp.dat <- tmp.dat[[1]]
      }
    }else{
      tmp.dat <- tmp.dat[[1]]
    }

    fcs.file.path <- file.path(sub("04_debarcoded", "05_concatenated", unique(dirname(fcs.paths.tmp))))

    if(!all(sapply(seq(length(input.cytofset)), function(i) length(input.cytofset[[i]]@description)))){
      "Missing keywords...?"
    }

    message("Building concatenate...")
    fcs.concatenated <- flowCore::fsApply(input.cytofset, flowCore::exprs)#concatenated expression matrix

    fcs.df <- data.frame(name = dimnames(fcs.concatenated)[[2]],
                         desc = tmp.dat$desc)
    fcs.df$range <- apply(apply(fcs.concatenated,2,range),2,diff)
    fcs.df$minRange <- apply(fcs.concatenated,2,min)
    fcs.df$maxRange <- apply(fcs.concatenated,2,max)

    df.anno <- Biobase::AnnotatedDataFrame(fcs.df)
    if(any(grepl("N", rownames(df.anno)))){
      rownames(df.anno) <- gsub("N", "", rownames(df.anno))
    }#need to sub out N or flowCore:::cols_to_pd will fail @ 'new_pid <- max(as.integer(gsub("\\$P", "", rownames(pd)))) + 1'

    fcs.concatenated <- methods::new("flowFrame",
                                     exprs = fcs.concatenated,
                                     parameters = df.anno,
                                     description = as.list(keyword.headers)
    )

    fcs.concatenated@description$POOL.ALIQUOT <- "CONCATENATED"

    aliquot.col <-  matrix(rep(seq(length(input.cytofset)),
                               flowCore::fsApply(input.cytofset, nrow)),
                           ncol = 1,
                           dimnames = list(NULL, c("aliquot_fcs")))

    fcs.concatenated <- flowCore::fr_append_cols(fcs.concatenated,aliquot.col)

    if(fcs.concatenated@description$`GLOBAL ID` == "UNASSIGNED"){
      fcs.concatenated@description$`$FIL` <- paste0(paste(fcs.concatenated@description$STUDY,
                                                          fcs.concatenated@description$POOL,
                                                          fcs.concatenated@description$CONDITION,
                                                          "UNASSIGNED",
                                                          "CONCATENATED",
                                                          sep = "_"),
                                                    ".fcs")
    }else if(fcs.concatenated@description$Visit == "Adult"){
      fcs.concatenated@description$`$FIL` <- paste0(paste(fcs.concatenated@description$Participant.Id,
                                                          fcs.concatenated@description$Visit,
                                                          fcs.concatenated@description$CONDITION,
                                                          fcs.concatenated@description$EXP,
                                                          "CONCATENATED",
                                                          sep = "_"),
                                                    ".fcs")
    }else{
      fcs.concatenated@description$`$FIL` <- paste0(paste(fcs.concatenated@description$Participant.Id,
                                                          paste0("V", fcs.concatenated@description$Sequence.Num),
                                                          fcs.concatenated@description$CONDITION,
                                                          "CONCATENATED",
                                                          sep = "_"),
                                                    ".fcs")
    }

    message(paste("Writing .fcs:", fcs.concatenated@description$`$FIL`, sep = " "))
    flowCore::write.FCS(fcs.concatenated, filename = file.path(fcs.file.path, fcs.concatenated@description$`$FIL`))

  }))
}

##asinh transform compensated .fcs expression data;default cofactor of 1 (divide by 2 to scale to fluor) for scatter; default cofactor of 500 for fluors
asinh.transform.fluor<-function(fcs,modified.fluor.cofactor=NULL,transform.scatter=F){
  if(transform.scatter){
    parms.scatter<-grep("SC",flowCore::colnames(fcs),value = T)
    fcs@exprs[,parms.scatter] <- asinh(fcs@exprs[,parms.scatter]/1)/2
  }
  m<-flowCore::markernames(fcs)
  parms.fluor<-names(m)
  if(!is.null(modified.fluor.cofactor)){
    m.mod <- names(modified.fluor.cofactor)
    parms.fluor.mod<-parms.fluor[!parms.fluor%in%names(grep(paste0(m.mod,collapse = "|"),m,value = T))]
    fcs@exprs[,parms.fluor.mod] <- asinh(fcs@exprs[,parms.fluor.mod]/500)
    for(i in m.mod){
      fcs@exprs[,names(grep(i,m,value = T))] <- asinh(fcs@exprs[,names(grep(i,m,value = T))]/modified.fluor.cofactor[[i]])
    }
  }else{
    fcs@exprs[,parms.fluor] <- asinh(fcs@exprs[,parms.fluor]/500)
  }
  return(fcs)
}
