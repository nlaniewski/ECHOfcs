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
