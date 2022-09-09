barcode.assignment.fsom.codes <- function(nCk = c(7,3), fsom.codes = NULL, threshold.override = NULL){

  if(nCk[1] != ncol(fsom.codes)){
    print("Wrong nCHOOSEk scheme!")
  }else{
    barcode.combinations <- t(utils::combn(nCk[1], nCk[2]))
    barcode.key <- matrix(data = 0,
                          nrow = nrow(barcode.combinations),
                          ncol = ncol(fsom.codes),
                          dimnames = list(NULL, colnames(fsom.codes)))

    for(i in seq(nrow(barcode.key))){
      barcode.key[i, barcode.combinations[i, ]] <- 1
    }
  }

  node.thresholds <- apply(fsom.codes, 2, function(x, high.edge.trim = "YES"){
    if(high.edge.trim == "YES"){
      q.val <- stats::quantile(x, probs = .999)
      d <- stats::density(x[x < q.val])
    }else{
      d <- stats::density(x)
    }
    valley.max.x = max(d$x[which(diff(sign(diff(d$y)))==2)])
    return(valley.max.x)
  })

  if(!is.null(threshold.override)){
    if(!is.list(threshold.override)){
      stop("Need a named list for override; example: 'list(Cd106Di = 0.75)'")
    }else{
      message(paste("Overring node threshold(s) for:", paste0(names(threshold.override), collapse = ","), sep = " "))
      for(i in names(threshold.override)){
        node.thresholds[[i]] <- threshold.override[[i]]
      }
    }
  }

  node.key <- sapply(names(node.thresholds), function(x){
    key <- ifelse(fsom.codes[, x] > node.thresholds[x], 1, 0)
  })

  node.key.t <- t(node.key)
  barcode.key.t <- t(barcode.key)

  barcode.assignment <- sapply(1:ncol(node.key.t), function(i) {
    if(max(colSums(node.key.t[, i] == barcode.key.t)) == nCk[1]){
      barcode <- which(colSums(node.key.t[, i] == barcode.key.t) == nCk[1])
    }else{
      barcode <- 0
    }
  })

  barcode.assignment <- list(threshold.type = "valley",
                             thresholds = node.thresholds,
                             barcodes = matrix(barcode.assignment, ncol = 1, dimnames = list(NULL, c("barcode"))),
                             key = barcode.key,
                             threshold.override = threshold.override
  )
  return(barcode.assignment)
}

mean.sd.bc.cut <- function(x){
  if(which.max(graphics::hist(x, plot = F)$counts)==1){
    bc.type <- 'bc.neg'
    sd.adjust <- 2.5
  }else{
    bc.type <- 'bc.pos'
    sd.adjust <- 3.5
  }

  dat.vec <- x[x>0]
  mean.vec <- mean(dat.vec)
  sd.vec <- stats::sd(dat.vec)

  cut <- c(mean.vec-(sd.vec*sd.adjust),
           mean.vec+(sd.vec*sd.adjust)
  )

  if(bc.type=='bc.neg'){
    cut <- cut[2]
  }
  cut <- cut[cut>0]
  if(length(cut)==1){
    cut.index <- which(x>cut[1])
  }else if(length(cut)==2){
    cut.index <- c(which(x<cut[1]),
                   which(x>cut[2]))
  }
  return(list(cut = cut,
              cut.index = cut.index))
}
##models a normal distribution based on vector input; used to generate cut values for the 'tails' in the original distribution
# qdist <- function(dat.vec, probs){
#   d <- density(dat.vec)
#   sd.dat.vec <- sd(dat.vec)
#   set.seed(20040501)
#   unname(quantile(rnorm(d$n,
#                         m = d$x[which.max(d$y)],
#                         sd = sd.dat.vec),
#                   probs = probs
#   ))
# }
##function to determine barcode cut 'bounds'; models a normal distribution then creates an index of events that fall outside the 99.99 percentile
# qdist.remove <- function(dat.bc){
#
#   barcode.val <- as.numeric(attributes(dat.bc)$comment)
#   barcode.positive <- combn(max(which(!colnames(dat.bc) %in% "barcode.col.index")),3)[,barcode.val]
#   barcode.negative <- which(!seq(1:max(which(!colnames(dat.bc) %in% "barcode.col.index"))) %in% combn(max(which(!colnames(dat.bc) %in% "barcode.col.index")),3)[,barcode.val])
#   barcode.names <- colnames(dat.bc)[!colnames(dat.bc) %in% "barcode.col.index"]
#   ##barcode.positive.bounds.remove
#   ##
#   barcode.positive.bounds <- setNames(vector(mode = "list", length = 3),
#                                       nm = barcode.names[barcode.positive])
#   for(i in names(barcode.positive.bounds)){
#     barcode.positive.bounds[[i]] <- qdist(dat.bc[[i]], probs = c(0.0001, 0.9999))#
#   }
#
#   barcode.positive.bounds.index <- setNames(vector(mode = "list", length = 3),
#                                             nm = names(barcode.positive.bounds))
#   for(i in names(barcode.positive.bounds.index)){
#     barcode.positive.bounds.index[[i]] <- c(which(dat.bc[[i]] < barcode.positive.bounds[[i]][1]),
#                                             which(dat.bc[[i]] > barcode.positive.bounds[[i]][2]))
#   }
#
#   barcode.positive.bounds.remove <- Reduce(union, barcode.positive.bounds.index)
#   ##barcode.negative.upper.bound.remove
#   ##
#   barcode.negative.upper.bound <- setNames(vector(mode = "numeric", length = max(which(!colnames(dat.bc) %in% "barcode.col.index"))-3),
#                                            nm = barcode.names[barcode.negative])
#   for(i in names(barcode.negative.upper.bound)){
#     barcode.negative.upper.bound[[i]] <- qdist(dat.bc[[i]][dat.bc[[i]] > 0.1], probs = 0.99)
#   }
#
#   barcode.negative.upper.bound.index <- setNames(vector(mode = "list", length = max(which(!colnames(dat.bc) %in% "barcode.col.index"))-3),
#                                                  nm = names(barcode.negative.upper.bound))
#   for(i in names(barcode.negative.upper.bound.index)){
#     barcode.negative.upper.bound.index[[i]] <- which(dat.bc[[i]] > barcode.negative.upper.bound[[i]])
#   }
#
#   barcode.negative.upper.bound.remove <- Reduce(union, barcode.negative.upper.bound.index)
#   ##qdist.remove
#   qdist.remove.index <- union(barcode.negative.upper.bound.remove, barcode.positive.bounds.remove)
#   qdist.remove.list <- list(index = union(barcode.negative.upper.bound.remove, barcode.positive.bounds.remove),
#                             bc.positive = barcode.positive.bounds,
#                             bc.negative = barcode.negative.upper.bound)
#   qdist.remove.list
# }
pair.names <- function(character.vec){
  if(!(length(character.vec)/2)%%1==0){
    pairs <- ceiling(length(character.vec)/2)
    pairs.max <- pairs*2
  }else if((length(character.vec)/2)%%1==0){
    pairs <- length(character.vec)/2
    pairs.max <- length(character.vec)
  }

  pair.1 <- seq(1, pairs.max-1, 2)
  pair.2 <- seq(2, pairs.max, 2)

  if(!(length(character.vec)/2)%%1==0){
    if(!all(pair.2 %in% seq(2, pairs.max-1, 2))){
      pair.2[which(!pair.2 %in% seq(2, pairs.max-1, 2))] <- pair.2[which(!pair.2 %in% seq(2, pairs.max-1, 2))-1]
    }
  }

  pair.names <- lapply(seq(pairs), function(j){
    c(character.vec[pair.1[j]], character.vec[pair.2[j]])
  })
}

barcode_CD45_fsom_ECHO_batch_condition <- function(echo.fcs.trimmed.paths, data.only = FALSE, parallel_SOM = TRUE,
                                                   generate.new = FALSE, cut.method = c('qdist','SD'), density.plots = FALSE, ...) {

  if(!is.vector(echo.fcs.trimmed.paths) & !is.character(echo.fcs.trimmed.paths)){
    stop("Need a vector of trimmed .fcs files (single-batch")
  }
  ##
  stim.conditions <- basename(unique(dirname(echo.fcs.trimmed.paths)))
  batch.date <- unique(stringr::str_extract(basename(echo.fcs.trimmed.paths),"[0-9]{8}"))
  if(length(stim.conditions)==2&length(batch.date)==1){
    print(paste("Generating barcodes for",paste("ECHO",batch.date,sep = "_"),":",stim.conditions))
  }else{
    stop("Need a unique ECHO batch ('yyyymmdd') with SEB and UNSTIM conditions...")
  }
  ##
  echo.fcs.trimmed.paths.condition <- sapply(stim.conditions,function(i) grep(i,echo.fcs.trimmed.paths,value = T))
  ##
  invisible(
    sapply(names(echo.fcs.trimmed.paths.condition),function(i){
      if(length(grep("fsom_coded",
                     list.files(sub("02_trimmed", "03_coded", unique(dirname(echo.fcs.trimmed.paths.condition[[i]]))),
                                pattern = ".rds")
      )
      )>=1 & generate.new == FALSE){
        stop(paste("fsom_coded already exists for",i,"..."))
      }else if(length(grep("fsom_coded",
                           list.files(sub("02_trimmed", "03_coded", unique(dirname(echo.fcs.trimmed.paths.condition[[i]]))),
                                      pattern = ".rds")
      )
      )>=1 & generate.new == TRUE){
        message(paste("fsom_coded already exists for",i,"generating a new file..."))
      }else{
        message(paste("creating a fsom_coded file for",i,"..."))
      }
    }
    )
  )
  ##
  for(i in names(echo.fcs.trimmed.paths.condition)){
    #sapply(names(echo.fcs.trimmed.paths.condition),function(i, cut.method = match.arg(cut.method)){
    fsom <- list()
    class(fsom) <- "FlowSOM"

    fcs.files <- echo.fcs.trimmed.paths.condition[[i]]

    rolling.row.index <- function(row.totals.vec){
      row.end <- cumsum(row.totals.vec)
      row.start <- (row.end - row.totals.vec)+1
      rolling.row <- matrix(c(row.start, row.end), ncol = 2, dimnames = list(names(row.end), c("start", "end")))
      rolling.row
    }

    fsom$metaData <- rolling.row.index(sapply(fcs.files, function(i) as.numeric(flowCore::read.FCSheader(i, keyword = "$TOT")), USE.NAMES = T))

    if(length(unique(dirname(fcs.files))) == 1){
      dir.mod <- sub("02_trimmed", "03_coded", unique(dirname(fcs.files)))
    }else{
      stop("is this a list of files from a single ECHO batch/condition?")
    }

    if(!dir.exists(dir.mod)){
      dir.create(dir.mod, recursive = T)
    }

    dir.mod.p <- file.path(dirname(dir.mod), "coded_plots")

    if(!dir.exists(dir.mod.p)){
      dir.create(dir.mod.p, recursive = T)
    }

    message(paste("Parallel reading of fcs files:",i))
    ##setup clusters for parallel computing
    cl <- parallel::makeCluster(parallel::detectCores())
    #cut.method <- match.arg(cut.method)
    if(cut.method == 'qdist'){
      parallel::clusterExport(cl, c("qdist", "qdist.remove"))
    }
    ##
    fcs.list <- parallel::parSapply(cl,fcs.files,function(fcs.file){
      ECHOfcs::read.fcs.selected.markers(fcs.path = fcs.file,
                                         selected.markers = "CD45",
                                         suppress.check = T)
    })
    parallel::stopCluster(cl);gc()
    fcs.list <- lapply(fcs.list,flowCore::exprs)
    ##
    fsom$data <- do.call(rbind,fcs.list);rm(fcs.list);gc()
    ##
    if(all(grepl("Cd", colnames(fsom$data)))){
      #generate Cadmium spill matrix
      cadmium.spill.values <- list(
        Cd106Di = stats::setNames(c(0.0, 0.2, 0.2, 0.6, 0.8, 0.4, 0.0),
                                  nm =paste0("Cd", c(106,110:114,116), "Di")
        ),
        Cd110Di = stats::setNames(c(0.0, 0.0, 1.4, 1.3, 0.5, 0.8, 0.2),
                                  nm =paste0("Cd", c(106,110:114,116), "Di")
        ),
        Cd111Di = stats::setNames(c(0.0, 0.3, 0.0, 4.5, 0.5, 0.8, 0.1),
                                  nm =paste0("Cd", c(106,110:114,116), "Di")
        ),
        Cd112Di = stats::setNames(c(0.0, 0.0, 0.9, 0.0, 1.0, 0.0, 0.0),
                                  nm =paste0("Cd", c(106,110:114,116), "Di")
        ),
        Cd113Di = stats::setNames(c(0.0, 0.1, 0.1, 1.5, 0.0, 4.0, 0.2),
                                  nm =paste0("Cd", c(106,110:114,116), "Di")
        ),
        Cd114Di = stats::setNames(c(0.0, 0.1, 0.1, 0.3, 0.3, 0.0, 0.2),
                                  nm =paste0("Cd", c(106,110:114,116), "Di")
        ),
        Cd116Di = stats::setNames(c(0.0, 0.3, 0.3, 0.7, 0.4, 1.8, 0.0),
                                  nm =paste0("Cd", c(106,110:114,116), "Di")
        )
      )
      A <- do.call(rbind, cadmium.spill.values)/100; diag(A) <- 1
      A <- t(A)
    }
    ##
    cl <- parallel::makeCluster(parallel::detectCores())
    parallel::clusterExport(cl, "A", envir = environment())
    ##
    message(paste("Parallel cadmium compensation of",i))
    ##
    fsom$data <- parallel::parRapply(cl, fsom$data, function(x){stats::coef(nnls::nnls(A = A, b = x))})
    parallel::stopCluster(cl);gc()
    ##
    fsom$data <- matrix(fsom$data,
                        nrow = length(fsom$data)/length(colnames(A)),
                        byrow = TRUE,
                        dimnames = list(NULL, colnames(A))
    )
    ##
    fsom$data <- apply(fsom$data, 2, function(x, cofactor = 5){ asinh(x)/cofactor})
    ##
    fsom$seed <- 20040501
    ##
    if(data.only == TRUE){
      return(fsom)
      stop("fsom.object with data only; no 'BuildSOM'")
    }
    ##
    if(parallel_SOM == TRUE){
      message(paste("Parallel SOM:",paste(ncol(fsom$data), "dimensions", sep = " "), sep = " "))
      message(paste0(colnames(fsom$data), collapse = " "))
      set.seed(fsom$seed)
      fsom$map <- EmbedSOM::SOM(fsom$data, xdim = 14, ydim = 14, batch=T, parallel=T, rlen=100)
    }else{
      message(paste("FlowSOM:",paste(ncol(fsom$data), "dimensions", sep = " "), sep = " "))
      message(paste0(colnames(fsom$data), collapse = " "))
      set.seed(fsom$seed)
      fsom <- FlowSOM::BuildSOM(fsom, xdim = 12, ydim = 12)
    }
    ##
    nCk <- c(ncol(fsom$data), 3)
    message(paste("Assigning barcodes based on a",paste(paste(nCk[1], "choose", "3", sep = "-"), "scheme", sep = " "), sep = " "))
    fsom$barcode.assignment <- barcode.assignment.fsom.codes(fsom.codes = fsom$map$codes, nCk = nCk, ...)
    ##
    if(cut.method == 'SD'){
      message("Cutting barcodes using 'SD' method...")

      b.col <- fsom$barcode.assignment$barcodes[fsom$map$mapping[, 1]]

      cut.list <- lapply(split(data.frame(fsom$data),b.col), function(i){
        if(any(grepl("threshold.override",names(fsom$barcode.assignment)))&length(fsom$barcode.assignment$threshold.override)>0){
          #message("threshold override")
          cut.list <- lapply(i[,!colnames(i) %in% names(fsom$barcode.assignment$threshold.override)], function(j){
            mean.sd.bc.cut(j)
          })
          if(length(fsom$barcode.assignment$threshold.override)==1){
            cut.name <- names(fsom$barcode.assignment$threshold.override)
            cut.list[[cut.name]]$cut <- fsom$barcode.assignment$threshold.override[[cut.name]]
            if(length(cut.list[[cut.name]]$cut)==1){
              cut.list[[cut.name]]$cut.index <- which(i[,cut.name]>fsom$barcode.assignment$threshold.override[[cut.name]])
            }else{
              stop("need to write code to handle two cuts (high/low); threshold override")
            }
          }else{
            stop("need to write code to handle threshold override of more than one channel")
          }
        }else{
          cut.list <- lapply(i, FUN = mean.sd.bc.cut)
        }
        cut.list$cut.index.union <- Reduce(union, lapply(cut.list, '[[', 2))
        return(cut.list)
      })
      cut.list <- cut.list[names(cut.list)!=0]

      barcodes <- as.numeric(names(cut.list))

      cut.zeroes <- unlist(sapply(barcodes, function(i){
        which(b.col==i)[cut.list[[as.name(i)]]$cut.index.union]
      }), use.names = F)

      fsom$barcode.col <- b.col
      fsom$barcode.col.cut <- fsom$barcode.col
      fsom$barcode.col.cut[cut.zeroes] <- 0
      fsom$barcodes.unique <- barcodes
      ##
      message("Generating barcode cut plots ('SD' method)...")
      plot.pairs <- pair.names(colnames(fsom$data))
      f.name <- unique(stringr::str_extract(rownames(fsom$metaData), "ECHO_[0-9]{8}_[A-Z]*"))
      ##
      cut.plots <- sapply(barcodes, function(i){
        dat.sub <- data.frame(fsom$data[fsom$barcode.col==i,])
        cuts <- sapply(cut.list[[as.name(i)]], '[', 'cut')
        bc.pos <- sub("Di.cut", "",names(which(sapply(cuts[!is.na(cuts)], length)==2)))
        plot.list <- sapply(plot.pairs, function(j){
          ggplot2::ggplot(dat.sub, ggplot2::aes(!!ggplot2::sym(j[1]),!!ggplot2::sym(j[2]))) +
            ggplot2::geom_hex(bins = 200) +
            viridis::scale_fill_viridis(option = "plasma", limits = c(0, 100), oob = scales::squish) +
            ggplot2::geom_vline(xintercept = cut.list[[as.name(i)]][[j[1]]]$cut) +
            ggplot2::geom_hline(yintercept = cut.list[[as.name(i)]][[j[2]]]$cut) +
            ggplot2::guides(fill = 'none') +
            ggplot2::xlim(-0.01,1.75) +
            ggplot2::ylim(-0.01,1.75)
        }, simplify = F)
        plot.page <- ggpubr::ggarrange(plotlist = plot.list)
        plot.page <- ggpubr::annotate_figure(plot.page,
                                             top = ggpubr::text_grob(paste(f.name,
                                                                           paste(paste("Barcode:", i, sep = " "), paste0(bc.pos, collapse = ",")),
                                                                           sep = "\n"),
                                                                     color = "black", face = "bold", size = 12)
        )
      }, simplify = F)
      ##
      out.path <- file.path(dir.mod.p,paste0(paste(f.name,"barcode_cuts_SD", Sys.Date(), sep = "_"),".pdf"))
      grDevices::pdf(out.path)
      print(cut.plots)
      grDevices::dev.off()
      #Biobase::openPDF(normalizePath(out.path))
    }else if(cut.method == 'qdist'){
      # message("Cutting barcodes using 'qdist' method...")
      # fsom$barcode.col <- fsom$barcode.assignment$barcodes[fsom$map$mapping[, 1]]
      # fsom$barcodes.unique <- setNames(as.integer(sort(unique(fsom$barcode.col[fsom$barcode.col>0]))),
      #                                  nm = paste0("barcode_", as.integer(sort(unique(fsom$barcode.col[fsom$barcode.col>0])))))
      #
      # fsom$barcode.cuts <- parLapply(cl, sapply(fsom$barcodes.unique, function(i) generate.barcoded.data.frame(fsom, barcode.val = i), simplify = F), function(i) qdist.remove(i))
      #
      # stopCluster(cl)
      #
      # new.zeroes <- unlist(sapply(names(fsom$barcodes.unique), function(i){
      #   which(fsom$barcode.col == fsom$barcodes.unique[[i]])[fsom$barcode.cuts[[i]]$index]
      # }), use.names = F)
      #
      # fsom$barcode.col.cut <- fsom$barcode.col
      # fsom$barcode.col.cut[new.zeroes] <- 0
      # ##
      # ##generate some plots
      # dat.bcs.plot.list <- sapply(names(fsom$barcodes.unique), function(i){
      #   dat.bc <- generate.barcoded.data.frame(fsom, barcode.val = fsom$barcodes.unique[[i]])
      #   bounds <- c(fsom$barcode.cuts[[i]]$bc.positive, fsom$barcode.cuts[[i]]$bc.negative)
      #   p.names <- pair.names(grep("barcode", colnames(dat.bc), invert = T, value = T))
      #
      #   plot.bc.pairs <- lapply(p.names, function(j){
      #     require(ggplot2)
      #     ggplot(dat.bc, aes_string(x = j[1], y = j[2])) +
      #       geom_hex(bins = 200, show.legend = F) +
      #       viridis::scale_fill_viridis(option = "plasma", limits = c(0, 100), oob = scales::squish) +
      #       xlim(-0.01,1.75) +
      #       ylim(-0.01,1.75) +
      #       geom_vline(xintercept = bounds[[j[1]]]) +
      #       geom_hline(yintercept = bounds[[j[2]]])
      #   })
      #
      #   attributes(plot.bc.pairs)$comment <- comment(dat.bc)
      #
      #   plot.bc.pairs
      # }, simplify = F)
      #
      # message("Generating barcode cut plots ('qdist' method)...")
      #
      # plot.bc.grobs <- sapply(dat.bcs.plot.list, function(i){
      #   require(gridExtra)
      #   barcode.val <- as.numeric(comment(i))
      #   b.name <- paste("barcode", barcode.val, sep = "_")
      #   if(ceiling(ncol(fsom$data)/2) == 4){
      #     plot.bc.grobs <- arrangeGrob(i[[1]],
      #                                  i[[2]],
      #                                  i[[3]],
      #                                  i[[4]],
      #                                  ncol = 2,
      #                                  top = paste(paste0(b.name, ":"),
      #                                              paste(combn(7,3)[, barcode.val], collapse = " "), sep = " ")
      #     )
      #   }else if(ceiling(ncol(fsom$data)/2) == 3){
      #     plot.bc.grobs <- arrangeGrob(i[[1]],
      #                                  i[[2]],
      #                                  i[[3]],
      #                                  ncol = 2,
      #                                  top = paste(paste0(b.name, ":"),
      #                                              paste(combn(6,3)[, barcode.val], collapse = " "), sep = " ")
      #     )
      #   }
      # }, simplify = F)
      #
      # rm(dat.bcs.plot.list)
      #
      # pool.name <- unique(sapply(basename(rownames(fsom$metaData)), function(i){
      #   paste0(strsplit(i, "_")[[1]][1:3], collapse = "_")
      # }))
      #
      # cond <- unique(sapply(basename(rownames(fsom$metaData)), function(i){
      #   paste0(strsplit(i, "_")[[1]][3], collapse = "_")
      # }))
      #
      # pdf(file.path(dir.mod.p, paste0(paste(pool.name, "barcode_cuts", Sys.Date(), sep = "_"), ".pdf")))
      # invisible(sapply(plot.bc.grobs, plot))
      # dev.off()
    }
    ##
    message(paste("Generating barcode total counts:",i))
    pool.name <- unique(stringr::str_extract(rownames(fsom$metaData), "ECHO_[0-9]{8}_[A-Z]*"))
    barcode_totals <- data.frame(stats::setNames(plyr::count(fsom$barcode.col.cut), nm = c("barcode", "count")),
                                 proportion = plyr::count(fsom$barcode.col.cut)$freq/sum(plyr::count(fsom$barcode.col.cut)$freq)*100,
                                 pool = pool.name,
                                 condition = strsplit(pool.name, "_")[[1]][3],
                                 aliquot = "CONCATENATED",
                                 experiment = sub("_", "", unique(stringr::str_extract(rownames(fsom$metaData), "_[0-9]{3}"))),
                                 date = unique(stringr::str_extract(rownames(fsom$metaData), "[0-9]{8}"))
    )
    barcode_totals$barcode <- as.factor(barcode_totals$barcode)
    barcode_totals.file.path <- file.path(dir.mod, "counts")
    ##
    if(!dir.exists(barcode_totals.file.path)){
      dir.create(barcode_totals.file.path, recursive = T)
    }
    ##
    barcode.counts.fname <- paste0(paste(pool.name, "counts", Sys.Date(), sep = "_"), ".csv")
    utils::write.csv(barcode_totals, file = file.path(barcode_totals.file.path, barcode.counts.fname), row.names = F)
    ##
    message(paste("Generating barcode yield plots:",i))
    dat <- droplevels(subset(barcode_totals, barcode_totals$barcode != 0))
    dat.zero <- droplevels(subset(barcode_totals, barcode_totals$barcode == 0))
    barcode_yields <- list(proportion = ggplot2::ggplot(dat, ggplot2::aes(!!quote(barcode), !!quote(proportion))) +
                             ggplot2::geom_col() +
                             ggplot2::labs(title = barcode_totals$pool,
                                           x = "Barcode",
                                           y = "Proportion (%) of Pool",
                                           caption = paste(paste0(round(dat.zero$proportion, 2),"%"), "unassigned", sep = " ")),
                           count = ggplot2::ggplot(dat, ggplot2::aes(!!quote(barcode),!!quote(count))) +
                             ggplot2::geom_col() +
                             ggplot2::labs(title = barcode_totals$pool,
                                           x = "Barcode",
                                           y = "Cell Count in Pool",
                                           caption = paste(dat.zero$count, "unassigned", sep = " "))
    )
    ##
    grDevices::pdf(file.path(dir.mod.p,paste0(paste(pool.name,"barcode_yields",Sys.Date(),sep = "_"),".pdf")))
    print(barcode_yields$proportion)
    print(barcode_yields$count)
    grDevices::dev.off()
    ##
    message(paste("Generating barcode medians:",i))
    barcode.medians <- t(sapply(stats::setNames(nm = sort(unique(fsom$barcode.col.cut))), function(i){
      apply(fsom$data[fsom$barcode.col.cut == i, ], 2, stats::median)
    }))

    barcode.medians.file.path <- file.path(dir.mod, "barcode_medians")
    if(!dir.exists(barcode.medians.file.path)){
      dir.create(barcode.medians.file.path, recursive = T)
    }

    barcode.medians.fname <- paste0(paste(pool.name, "barcode_medians", Sys.Date(), sep = "_"), ".csv")
    utils::write.csv(barcode.medians, file = file.path(barcode.medians.file.path, barcode.medians.fname), row.names = T)

    pheatmap::pheatmap(barcode.medians,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       main = paste(pool.name,"\nBarcode Assignment: Median Intensities"),
                       filename = file.path(dir.mod.p, paste0(paste(pool.name, "barcode_medians", Sys.Date(), sep = "_"), ".pdf"))
    )
    ##
    if(density.plots){
      # message("Generating density plots...")
      # xlims <- c(0, max(apply(fsom$data[fsom$barcode.col.cut != 0, ], 2, function(x) quantile(x, probs = 0.9999))))
      # xlims[2] <- xlims[2]+xlims[2]*.005
      #
      # dat.densities <- lapply(setNames(nm = sort(unique(fsom$barcode.col.cut))), function(i){
      #   dat.melt <- reshape2::melt(as.data.frame(fsom$data[fsom$barcode.col.cut == i, ]),
      #                              value.name = "expression",
      #                              variable.name = "marker")
      #   if(i == 0){
      #     dat.melt$code <- factor(0, levels = c(0,1))
      #   }else{
      #     dat.melt$code <- as.factor(rep(fsom$barcode.assignment$key[i, ], each = nrow(dat.melt)/ncol(fsom$barcode.assignment$key)))
      #   }
      #
      #   ggplot(dat.melt, aes(x = expression, fill = code)) +
      #     geom_density() +
      #     facet_wrap(~marker, ncol = 1, scales = "free_y") +
      #     xlim(xlims) +
      #     labs(title = pool.name,
      #          subtitle = paste("Barcode:", i, sep = " "))
      #
      # })
      #
      # pdf(file.path(dir.mod.p, paste0(paste(pool.name, "barcode_densities", Sys.Date(), sep = "_"), ".pdf")), width = 5, height = 10)
      # invisible(print(dat.densities))
      # dev.off()
    }
    ##
    fsom.fname <- paste("ECHO",
                        unique(sapply(basename(fcs.files), function(i) strsplit(i, "_")[[1]][2])),
                        unique(sapply(basename(fcs.files), function(i) strsplit(i, "_")[[1]][3])),
                        "fsom",
                        "coded",
                        sep = "_"
    )

    message("Saving fsom...")
    saveRDS(fsom, file.path(dir.mod, paste0(paste(fsom.fname, Sys.Date(), sep = "_"), ".rds")))
    message("fsom saved")
    #})
  }
}

fsom_coded_to_fcs_coded <- function(fsom_coded.file.path){
  message(paste("Reading:", fsom_coded.file.path))
  fsom <- readRDS(fsom_coded.file.path)
  if(length(which((names(fsom) %in% c("barcode.col","barcode.col.cut"))))!=2){
    stop("Need fsom.object with '$barcode.col' and '$barcode.col.cut'")
  }
  invisible(sapply(rownames(fsom$metaData), function(i){

    message(paste("Reading file:", i, sep = " "))

    fcs <- flowCore::read.FCS(i, transformation = F, truncate_max_range = F)

    samp.index <- which(rownames(fsom$metaData) == i)
    row.index <- fsom$metaData[samp.index,1]:fsom$metaData[samp.index,2]

    if(nrow(fcs) != length(row.index)){
      stop("Row total mismatch...")
    }

    fcs <- flowCore::fr_append_cols(fcs,
                                    matrix(c(fsom$map$mapping[, 1][row.index],
                                             fsom$barcode.col[row.index],
                                             fsom$barcode.col.cut[row.index]),
                                           ncol = 3,
                                           dimnames = list(NULL, c("node_CD45",
                                                                   "barcode",
                                                                   "barcode_cut"))
                                    )
    )

    f.name <- sub(".fcs", "", fcs@description$`$FIL`)

    fcs@description$`$FIL` <- paste0(paste(f.name, "fsom_coded", sep = "_"), ".fcs")

    fcs.file.path <- file.path(sub("02_trimmed", "03_coded", dirname(rownames(fsom$metaData)[1])))

    if(!dir.exists(fcs.file.path)){
      dir.create(fcs.file.path)
    }

    flowCore::write.FCS(fcs, filename = file.path(fcs.file.path, fcs@description$`$FIL`))
  }))
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

debarcode_batched <- function(coded.fcs.file.path, use.batched.sheets = TRUE, metadata.sheet.path = NULL, use.barcode.cut = TRUE){

  keyword.headers <- fcs.headers.no.pars(fcs.file.paths = coded.fcs.file.path, return.combined = T)

  if(!any(names(keyword.headers) %in% "EXP")){
    stop("need/expect 'EXP' keyword")
  }else if(!any(names(keyword.headers) %in% "STUDY")){
    stop("need/expect 'STUDY' keyword")
  }else if(!any(names(keyword.headers) %in% "POOL")){
    stop("need/expect 'POOL' keyword")
  }else if(!any(names(keyword.headers) %in% "CONDITION")){
    stop("need/expect 'CONDITION' keyword")
  }else if(!any(names(keyword.headers) %in% "POOL.ALIQUOT")){
    stop("need/expect 'POOL.ALIQUOT' keyword")
  }

  if(use.batched.sheets == TRUE){

    batched_merged.file.path <- grep(paste0(keyword.headers[c("STUDY", "EXP")], collapse = "_"),
                                     list.files("../metadata/004_batched_merged/", full.names = T, pattern = ".xlsx"), value = T)
    mdat <- readxl::read_xlsx(batched_merged.file.path)
    #NA fix for some of the earlier sheets...
    mdat$barcode.SEB_UNSTIM <- sub(",NA", "", mdat$barcode.SEB_UNSTIM)
    mdat$barcode.seb <- sapply(strsplit(mdat$barcode.SEB_UNSTIM, ","), '[', 1)
    mdat$barcode.unstim <- sapply(strsplit(mdat$barcode.SEB_UNSTIM, ","), '[', 2)
    mdat$name.seb <- paste("ECHO",unique(mdat$pool.date),"SEB",mdat$barcode.seb, sep = "_")
    mdat$name.unstim <- paste("ECHO",unique(mdat$pool.date),"UNSTIM",mdat$barcode.unstim, sep = "_")

    if(!all(c("GID.1", "Visit", "Participant.Id", "Sequence.Num") %in% colnames(mdat))){
      stop("Check mdat sheet; need 'GID', 'Visit', 'Participant.Id', 'Sequence.Num'...")
    }

    if(keyword.headers["POOL"] != unique(mdat$pool.date)){
      stop("Pool date mismatch...")
    }

    if(keyword.headers["CONDITION"] == "SEB"){
      named.barcodes <- stats::setNames(sapply(mdat$barcode.SEB_UNSTIM, function(i) as.numeric(strsplit(i, ",")[[1]][1]), USE.NAMES = F),
                                        nm = mdat$participant_sequence)
    }else if(keyword.headers["CONDITION"] == "UNSTIM"){
      named.barcodes <- stats::setNames(sapply(mdat$barcode.SEB_UNSTIM, function(i) as.numeric(strsplit(i, ",")[[1]][2]), USE.NAMES = F),
                                        nm = mdat$participant_sequence)
    }
    named.barcodes <- named.barcodes[!is.na(named.barcodes)]
    names(named.barcodes)[grep("HD", names(named.barcodes))] <- paste(grep("HD", names(named.barcodes), value = T),
                                                                      unique(mdat$pool.date),
                                                                      sep = "_")

  }else if(use.batched.sheets == FALSE & !is.null(metadata.sheet.path)){
    mdat <- readxl::read_excel(metadata.sheet.path, sheet = grep("manifest_batched", readxl::excel_sheets(metadata.sheet.path)))
    mdat <- mdat[which(mdat$batch.number == as.numeric(keyword.headers["EXP"])), ]
    if(nrow(mdat) == 0){
      stop("No batch found in metadata sheet...")
    }
  }else{
    stop("Need metadata sheet for assigning subject/global ids...")
  }

  dir.mod <- sub("03_coded", "04_debarcoded", unique(dirname(coded.fcs.file.path)))

  invisible(
    sapply(dir.mod, function(i){
      if(!dir.exists(i)){
        dir.create(i, recursive = T)
      }
    })
  )
  ##
  message(paste("Debarcoding:", coded.fcs.file.path, sep = " "))
  fcs <- flowCore::read.FCS(coded.fcs.file.path, transformation = F, truncate_max_range = F)

  if(use.barcode.cut == TRUE){
    barcode.vec <- "barcode_cut"
  }else{
    barcode.vec <- "barcode"
  }

  sapply(sort(unique(fcs@exprs[, barcode.vec])), function(i){

    fcs.tmp <- fcs[fcs@exprs[, barcode.vec] == i, ]

    fcs.tmp@description$BARCODE <- i

    if(fcs.tmp@description$BARCODE == 0){

      fcs.merge.name <- paste(fcs.tmp@description$STUDY,
                              fcs.tmp@description$POOL,
                              fcs.tmp@description$CONDITION,
                              "UNASSIGNED",
                              fcs.tmp@description$POOL.ALIQUOT,
                              sep = "_")

      fcs.tmp@description$`$FIL` <- paste0(fcs.merge.name,".fcs")

      fcs.tmp@description$`Participant Id` <- "UNASSIGNED"
      fcs.tmp@description$`GLOBAL ID` <- "UNASSIGNED"
      fcs.tmp@description$Visit <- "UNASSIGNED"
      fcs.tmp@description$`Sequence Num` <- "UNASSIGNED"
      fcs.tmp@description$date.thawed <- unique(mdat$date.thawed)

      flowCore::write.FCS(fcs.tmp, filename = file.path(dir.mod, fcs.tmp@description$`$FIL`))

    }else if(fcs.tmp@description$CONDITION == "UNSTIM"&!fcs.tmp@description$BARCODE %in% as.numeric(mdat$barcode.unstim[!is.na(mdat$barcode.unstim)])){
      message(paste("Issue with barcode:", i, sep = " "))
      message("erroneous 'UNSTIM' barcode/name; assigning these events to 'mismatch' file")
      fcs.merge.name <- paste(fcs.tmp@description$STUDY,
                              fcs.tmp@description$POOL,
                              fcs.tmp@description$CONDITION,
                              paste("MISMATCH", fcs.tmp@description$BARCODE, sep = "_"),
                              fcs.tmp@description$POOL.ALIQUOT,
                              sep = "_")

      fcs.tmp@description$`$FIL` <- paste0(fcs.merge.name,".fcs")

      fcs.tmp@description$`Participant Id` <- "MISMATCH"
      fcs.tmp@description$`GLOBAL ID` <- "MISMATCH"
      fcs.tmp@description$Visit <- "MISMATCH"
      fcs.tmp@description$`Sequence Num` <- "MISMATCH"
      fcs.tmp@description$date.thawed <- unique(mdat$date.thawed)

      flowCore::write.FCS(fcs.tmp, filename = file.path(dir.mod, fcs.tmp@description$`$FIL`))

    }else if(fcs.tmp@description$CONDITION == "SEB"&!fcs.tmp@description$BARCODE %in% as.numeric(mdat$barcode.seb[!is.na(mdat$barcode.seb)])){
      message(paste("Issue with barcode:", i, sep = " "))
      message("erroneous 'SEB' barcode/name; assigning these events to 'mismatch' file")
      fcs.merge.name <- paste(fcs.tmp@description$STUDY,
                              fcs.tmp@description$POOL,
                              fcs.tmp@description$CONDITION,
                              paste("MISMATCH", fcs.tmp@description$BARCODE, sep = "_"),
                              fcs.tmp@description$POOL.ALIQUOT,
                              sep = "_")

      fcs.tmp@description$`$FIL` <- paste0(fcs.merge.name,".fcs")

      fcs.tmp@description$`Participant Id` <- "MISMATCH"
      fcs.tmp@description$`GLOBAL ID` <- "MISMATCH"
      fcs.tmp@description$Visit <- "MISMATCH"
      fcs.tmp@description$`Sequence Num` <- "MISMATCH"
      fcs.tmp@description$date.thawed <- unique(mdat$date.thawed)

      flowCore::write.FCS(fcs.tmp, filename = file.path(dir.mod, fcs.tmp@description$`$FIL`))

    }else{
      fcs.merge.name <- paste(fcs.tmp@description$STUDY,
                              fcs.tmp@description$POOL,
                              fcs.tmp@description$CONDITION,
                              fcs.tmp@description$BARCODE,
                              sep = "_")

      if(fcs.tmp@description$CONDITION == "UNSTIM"){
        mdat.tmp <- mdat[which(mdat$name.unstim %in% fcs.merge.name), ]
        mdat.tmp <- mdat.tmp[, grep("seb", colnames(mdat.tmp), invert = T)]
      }else if(fcs.tmp@description$CONDITION == "SEB"){
        mdat.tmp <- mdat[which(mdat$name.seb %in% fcs.merge.name), ]
        mdat.tmp <- mdat.tmp[, grep("unstim", colnames(mdat.tmp), invert = T)]
      }

      for(j in seq(ncol(mdat.tmp))){
        fcs.tmp@description[colnames(mdat.tmp)[j]] <- mdat.tmp[[j]]
      }

      fcs.tmp@description$`GLOBAL ID` <- mdat.tmp$GID.1

      if(fcs.tmp@description$Visit == "Adult"){
        fcs.tmp@description$`$FIL` <- paste0(paste(fcs.tmp@description[[grep("participant.*id", colnames(mdat.tmp), value = T, ignore.case = T)]],
                                                   fcs.tmp@description[[grep("Visit", colnames(mdat.tmp), value = T, ignore.case = T)]],
                                                   fcs.tmp@description$CONDITION,
                                                   fcs.tmp@description$EXP,
                                                   fcs.tmp@description$POOL.ALIQUOT,
                                                   sep = "_"),
                                             ".fcs")
      }else{
        fcs.tmp@description$`$FIL` <- paste0(paste(fcs.tmp@description[[grep("participant.*id", colnames(mdat.tmp), value = T, ignore.case = T)]],
                                                   paste0("V", fcs.tmp@description[[grep("sequence.*num", colnames(mdat.tmp), value = T, ignore.case = T)]]),
                                                   fcs.tmp@description$CONDITION,
                                                   fcs.tmp@description$POOL.ALIQUOT,
                                                   sep = "_"),
                                             ".fcs")
      }

      flowCore::write.FCS(fcs.tmp, filename = file.path(dir.mod, fcs.tmp@description$`$FIL`))
    }
  })
}
