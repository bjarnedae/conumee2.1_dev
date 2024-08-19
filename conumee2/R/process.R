##### PROCESSING methods #####

#' CNV.fit
#' @description Normalize query sample intensities by fitting intensities to reference set using a linear regression model.
#' @param query \code{CNV.data} object of query sample (multiple samples).
#' @param ref \code{CNV.data} object of reference set.
#' @param anno \code{CNV.anno} object. Use \code{CNV.create_anno} do create.
#' @param intercept logical. Should intercept be considered? Defaults to \code{TRUE}.
#' @param ... Additional parameters (\code{CNV.fit} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The log2 ratio of query intensities versus a linear combination of reference set intensities that best reflects query intensities is calculated (as determined by linear regression). Every query sample in the \code{CNV.analysis} object is compared to the linear combination fo control samples individually. The annotations provided to \code{CNV.fit} are saved within the returned \code{CNV.analysis} object and used for subsequent analysis steps.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d, ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' #x <- CNV.bin(x)
#' #x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#'
#' # general information
#' x
#' show(x)
#'
#' # coefficients from linear regression
#' x@@fit$coef
#'
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.fit", function(query, ref, anno, ...) {
  standardGeneric("CNV.fit")
})

#' @rdname CNV.fit
setMethod("CNV.fit", signature(query = "CNV.data", ref = "CNV.data", anno = "CNV.anno"),
          function(query, ref, anno, intercept = TRUE) {
            if (length(query@intensity) == 0)
              stop("query intensities unavailable, run CNV.load")
            if (length(ref@intensity) == 0)
              stop("reference set intensities unavailable, run CNV.load")

            if (length(query@intensity) != 1)
              message("using multiple query samples")
            if (length(ref@intensity) == 1)
              warning("reference set contains only a single sample. use more samples for better results.")

            query.probes.pr <- as.logical() # check if all probes from query samples are within the annotation object
            for(i in 1:length(query@intensity)){
              pr <- !(all(is.element(names(anno@probes), names(query@intensity[[i]]))) | all(is.element(mcols(anno@probes)$EPICv2_Loci, names(query@intensity[[i]]))))
              query.probes.pr <- c(query.probes.pr, pr)
            }
            if(any(query.probes.pr)){
              stop(paste("query intensities not given for all probes, sample(s)", which(query.probes.pr), sep = " "))
            }

            ref.probes.pr <- as.logical() # check if all probes from reference samples are within the annotation object
            for(i in 1:length(ref@intensity)){
              pr <- !(all(is.element(names(anno@probes), names(ref@intensity[[i]]))) | all(is.element(mcols(anno@probes)$EPICv2_Loci, names(ref@intensity[[i]]))))
              ref.probes.pr <- c(ref.probes.pr, pr)
            }
            if(any(ref.probes.pr)){
              stop(paste("reference intensities not given for all probes, sample(s)", which(ref.probes.pr), sep = " "))
            }

            query.df <- CNV.df(query, anno)
            ref.df <- CNV.df(ref, anno)

            if(!all(substr(rownames(query.df), start = 1, stop = 10) == substr(rownames(ref.df), start = 1, stop = 10))){
              stop("probe ids for query and reference samples are not matching.")
            }

            object <- new("CNV.analysis")
            object@name <- colnames(query.df)
            object@date <- date()
            object@fit$args <- list(intercept = intercept)

            object@anno <- anno

            object@fit$coef <- data.frame(matrix(ncol = 0, nrow = ncol(ref.df)))
            object@fit$ratio <- data.frame(matrix(ncol = 0, nrow = nrow(query.df)))
            for (i in 1:ncol(query.df)) {

              message(paste(colnames(query.df)[i]), " (",round(i/ncol(query.df)*100, digits = 3), "%", ")", sep = "")
              r <- cor(query.df, ref.df)[i, ] < 0.99
              if (any(!r)) message("query sample seems to also be in the reference set. not used for fit.")
              if (intercept) {
                ref.fit <- lm(y ~ ., data = data.frame(y = log2(query.df[,i]), X = log2(ref.df[,r])))
              } else {
                ref.fit <- lm(y ~ . - 1, data = data.frame(y = log2(query.df[,i]), X = log2(ref.df[,r])))
              }
              coefs <- rep(NA, ncol(ref.df))
              coefs[r] <- as.numeric(ref.fit$coefficients[-1])
              object@fit$coef <- cbind(object@fit$coef,coefs)

              ref.predict <- predict(ref.fit)
              ref.predict[ref.predict < 0] <- 0

              object@fit$ratio <- cbind(object@fit$ratio, log2(query.df[,i]) - ref.predict)
            }

            colnames(object@fit$coef) <- colnames(query.df)
            rownames(object@fit$coef) <- colnames(ref.df)
            colnames(object@fit$ratio) <- colnames(query.df)
            rownames(object@fit$ratio) <- rownames(query.df)

            object@fit$noise <- as.numeric()
            for (i in 1:ncol(query.df)) {
              object@fit$noise <- c(object@fit$noise, sqrt(mean((object@fit$ratio[-1,i] - object@fit$ratio[-nrow(object@fit$ratio),i])^2,na.rm = TRUE)))
            }

            names(object@fit$noise) <- colnames(query.df)
            return(object)
          })

#' CNV.df
#' @description Turn the probe intensities within a \code{CNV.data} object (list) into a dataframe.
#' @param object \code{CNV.data} object.
#' @param anno \code{CNV.anno} object.
#' @return \code{CNV.data} object.
#' @author Bjarne Daenekas \email{conumee@@hovestadt.bio}
setGeneric("CNV.df", function(object, anno) {
  standardGeneric("CNV.df")
})

#' @rdname CNV.df
setMethod("CNV.df", signature(object = "CNV.data"), function(object, anno) {

  probes <- anno@probes
  if(all(nchar(lapply(object@intensity, function(x) names(x)[1])) == 15)){
    p <- mcols(probes)$EPICv2_Loci
    df.object <- as.data.frame(matrix(nrow = length(p) , ncol = 0))
    for(i in 1:length(object@intensity)){
      object.p <- object@intensity[[i]][p]
      df.object <- cbind(df.object, object.p)
    }
    colnames(df.object) <- names(object@intensity)
  } else{
    df.object <- as.data.frame(matrix(nrow = length(probes) , ncol = 0))
    for(i in 1:length(object@intensity)){
      if(nchar(names(object@intensity[[i]][1])) == 15){
        p <- mcols(probes)$EPICv2_Loci
        object.p <- object@intensity[[i]][p]
        names(object.p) <- substr(names(object.p), start = 1, stop = 10)
      }else{
        p <- names(probes)
        object.p <- object@intensity[[i]][p]
      }
      df.object <- cbind(df.object, object.p)
    }
    colnames(df.object) <- names(object@intensity)
  }
  return(df.object)
})

#' CNV.bin
#' @description Combine single probe intensitiy values into predefined bins.
#' @param object \code{CNV.analysis} object.
#' @param ... Additional parameters (\code{CNV.bin} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The median intensity per bin and its variance are calculated. Bins are defined using \code{CNV.create_anno}. A value by which all probe and bin intensity values are shifted in subsequent analysis steps is calculated by minimizing the median absolute deviation from all bins to zero (ideally shifting the copy-neutral state to 0).
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' x <- CNV.bin(x)
#' #x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#'
#' # general information
#' x
#' show(x)
#'
#' # coefficients of linear regression
#' coef(x)
#'
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.bin", function(object, ...) {
  standardGeneric("CNV.bin")
})

#' @rdname CNV.bin
setMethod("CNV.bin", signature(object = "CNV.analysis"), function(object) {
  if (length(object@fit) == 0)
    stop("fit unavailable, run CNV.fit")

  if(all(nchar(rownames(object@fit$ratio)) == 15)){
    n.probes <- mcols(object@anno@probes)$EPICv2_Loci
  }else{
    n.probes <- names(object@anno@probes)
  }

  o1 <- as.matrix(findOverlaps(query = object@anno@bins, subject = object@anno@probes))
  o2 <- data.frame(bin = names(object@anno@bins)[o1[, "queryHits"]],
                   probe = n.probes[o1[, "subjectHits"]], stringsAsFactors = FALSE)

  object@bin$ratio <- vector(mode = "list", length = ncol(object@fit$ratio))
  object@bin$variance <- vector(mode = "list", length = ncol(object@fit$ratio))
  object@bin$shift <- as.numeric()
  for (i in 1:ncol(object@fit$ratio)) {

    message(paste(colnames(object@fit$ratio)[i]), " (",round(i/ncol(object@fit$ratio)*100, digits = 3), "%", ")", sep = "")

    object@bin$ratio[[i]] <- sapply(split(object@fit$ratio[o2[, "probe"],i], o2[,"bin"]),
                                    median, na.rm = TRUE)[names(object@anno@bins)]

    object@bin$variance[[i]] <- sapply(split(object@fit$ratio[o2[, "probe"],i], o2[,"bin"]),
                                       var, na.rm = TRUE)[names(object@anno@bins)]

    object@bin$shift <- c(object@bin$shift, optim(0, function(s) median(abs(object@bin$ratio[[i]] -s),na.rm = TRUE),
                                                  method = "Brent", lower = -100, upper = 100)$par)
  }

  names(object@bin$shift) <- colnames(object@fit$ratio)
  names(object@bin$ratio) <- colnames(object@fit$ratio)
  names(object@bin$variance) <- colnames(object@fit$ratio)

  return(object)
})

#' CNV.detail
#' @description Combine single probe values within detail regions.
#' @param object \code{CNV.analysis} object.
#' @param ... Additional parameters (\code{CNV.detail} generic, currently not used).
#' @return \code{CNV.analysis} object.
#' @details The median intensity per detail region is calculated. Detail regions are defined using \code{CNV.create_anno(detail_bed=)}
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' x <- CNV.bin(x)
#' x <- CNV.detail(x)
#' #x <- CNV.segment(x)
#'
#' # general information
#' x
#' show(x)
#'
#' # coefficients of linear regression
#' coef(x)
#'
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.detail", function(object, ...) {
  standardGeneric("CNV.detail")
})

#' @rdname CNV.detail
setMethod("CNV.detail", signature(object = "CNV.analysis"), function(object) {
  if (length(object@fit) == 0)
    stop("fit unavailable, run CNV.fit")
  # if(length(object@bin) == 0) stop('bin unavailable, run CNV.bin')

  if (length(object@anno@detail) == 0) {
    message("no detail regions provided, define using CNV.create_anno")
  } else {

    if(all(nchar(rownames(object@fit$ratio)) == 15)){
      n.probes <- mcols(object@anno@probes)$EPICv2_Loci
    }else{
      n.probes <- names(object@anno@probes)
    }

    d1 <- as.matrix(findOverlaps(query = object@anno@detail, subject = object@anno@probes))
    d2 <- data.frame(detail = values(object@anno@detail)$name[d1[,"queryHits"]], probe = n.probes[d1[, "subjectHits"]], stringsAsFactors = FALSE)

    object@detail$ratio <- vector(mode = "list", length = ncol(object@fit$ratio))
    for (i in 1:ncol(object@fit$ratio)) {
      object@detail$ratio [[i]]<- sapply(split(object@fit$ratio[d2[, "probe"],i],
                                               d2[, "detail"]), median, na.rm = TRUE)[values(object@anno@detail)$name]
    }
    names(object@detail$ratio) <- colnames(object@fit$ratio)
    object@detail$n_probes <- table(d2[, 1])[values(object@anno@detail)$name]

  }
  return(object)
})

#' @import DNAcopy
NULL

#' CNV.segment
#' @description Segment bin values (wrapper of \code{DNAcopy} package). Each bin is assigned to a weight that is inversely proprotional to its variance.
#' @param object \code{CNV.analysis} object.
#' @param alpha See details. Defaults to 0.001.
#' @param nperm See details. Defaults to 50000.
#' @param min.width See details. Defaults to 5.
#' @param undo.splits See details. Defaults to 'sdundo'.
#' @param undo.SD See details. Defaults to 2.2.
#' @param verbose See details. Defaults to 0.
#' @param ... Additional parameters supplied to the \code{segment} method of the \code{DNAcopy} package.
#' @return \code{CNV.analysis} object.
#' @details This method is a wrapper of the CNA, segment, segments.summary and segments.p methods of the DNAcopy package. Please refer to the respective man pages for more detailed information. The default parameters of \code{CNV.segment} override some of the default parameters of segment and are optimized for 450k data CNV analysis.
#' @examples
#' # prepare
#' library(minfiData)
#' data(MsetEx)
#' d <- CNV.load(MsetEx)
#' data(detail_regions)
#' anno <- CNV.create_anno(detail_regions = detail_regions)
#'
#' # create object
#' x <- CNV.fit(query = d['GroupB_1'], ref = d[c('GroupA_1', 'GroupA_2', 'GroupA_3')], anno)
#'
#' # modify object
#' x <- CNV.bin(x)
#' x <- CNV.detail(x)
#' x <- CNV.segment(x)
#'
#' # general information
#' x
#' show(x)
#'
#' # coefficients of linear regression
#' coef(x)
#'
#' # show or replace sample name
#' names(x)
#' names(x) <- 'Sample 1'
#' @author Volker Hovestadt, Bjarne Daenekas \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.segment", function(object, ...) {
  standardGeneric("CNV.segment")
})

#' @rdname CNV.segment
setMethod("CNV.segment", signature(object = "CNV.analysis"), function(object,
                                                                      alpha = 0.001, nperm = 50000, min.width = 5, undo.splits = "sdundo",
                                                                      undo.SD = 2.2, verbose = 0, ...) {
  if(length(object@fit) == 0){
    stop('fit unavailable, run CNV.fit')
  }
  if (length(object@bin) == 0){
    stop("bin unavailable, run CNV.bin")
  }

  a1 <- formals()
  a2 <- as.list(match.call())[-1]
  object@seg$args <- as.list(sapply(setdiff(unique(names(c(a1, a2))),
                                            c("object", "verbose")), function(an) if (is.element(an, names(a2)))
                                              a2[[an]] else a1[[an]], simplify = FALSE))

  object@seg$summary <- vector(mode = "list", length = ncol(object@fit$ratio))
  object@seg$p <- vector(mode = "list", length = ncol(object@fit$ratio))

  for (i in 1:ncol(object@fit$ratio)) {

    message(paste(colnames(object@fit$ratio)[i]), " (",round(i/ncol(object@fit$ratio)*100, digits = 3), "%", ")", sep = "")

    x1 <- DNAcopy::CNA(genomdat = object@bin$ratio[[i]][names(object@anno@bins)],
                       chrom = as.vector(seqnames(object@anno@bins)), maploc = values(object@anno@bins)$midpoint,
                       data.type = "logratio", sampleid = "sampleid")

    x2 <- DNAcopy::segment(x = x1, weights = 1/object@bin$variance[[i]][names(object@anno@bins)], verbose = verbose, min.width = min.width,
                           nperm = nperm, alpha = alpha, undo.splits = undo.splits, undo.SD = undo.SD,
                           ...)

    object@seg$summary[[i]] <- DNAcopy::segments.summary(x2)
    object@seg$summary[[i]]$chrom <- as.vector(object@seg$summary[[i]]$chrom)
    object@seg$summary[[i]]$ID <- colnames(object@fit$ratio)[i]
    object@seg$p[[i]] <- DNAcopy::segments.p(x2)
    object@seg$p[[i]]$chrom <- as.vector(object@seg$p[[i]]$chrom)
  }
  names(object@seg$summary) <- colnames(object@fit$ratio)
  names(object@seg$p) <- colnames(object@fit$ratio)
  return(object)
})

#' @import nullranges
NULL

#' CNV.focal
#' @description This function aims to detect focal CNVs that are defined as \code{detail_regions} in the annotation object.
#' @param object \code{CNV.analysis} object.
#' @param sig_cgenes logical. Should the genes from the Cancer Gene Census be assessed as well? (The high number of genes can lead to false positive results.) Default to \code{FALSE}.
#' @param conf numeric. Confidence level to determine the log2-threshold for focal alterations. Default to \code{0.99}.
#' @param R numeric. Parameter for the \code{bootRanges} function. The number of bootstrap samples to generate. Default to \code{100}.
#' @param blockLength numeric. Parameter for the \code{bootRanges} function. The length (in basepairs) of the blocks for Segmented Block Bootstrapping. Default to \code{500000}.
#' @param proportionLength logical. From the \code{nullranges} package: for the segmented block bootstrap, whether to use scaled block lengths, (scaling by the proportion of the segmentation state out of the total genome length)
#' @param ... Additional parameters (\code{CNV.detailplot} generic, currently not used).
#' @return A \code{CNV.analysis} object with significantly altered bins and predefined focal regions of interest.
#' @details This function should facilitate the detection of CNVs that affect single genes. K-means clustering is used to assign each segment to a copy-number state. The mean and standard deviation for each state derived from Segmented Block Bootstrapping is used to model a normal distribution for each state. The confidence interval is used to define state-spefici thresholds for focal CNVs.
#' @examples
#'
#' x <- CNV.focal(x, sig_cgenes = TRUE, conf = 0.99, R = 100, blockLength = 500000)
#' x@@detail$del.bins  #bins that are part of deletions.
#' x@@detail$amp.bins  #bins that are part of amplifications.
#' x@@detail$del.detail.regions #deleted predefined regions
#' x@@detail$amp.detail.regions #amplified predefined regions
#' x@@detail$del.cancer.genes #deleted genes from the Cancer Gene Census
#' x@@detail$amp.cancer.genes #amplified genes from the Cancer Gene Census
#'
#'
#' @author Bjarne Daenekas, Volker Hovestadt \email{conumee@@hovestadt.bio}
#' @export
setGeneric("CNV.focal", function(object, ...) {
  standardGeneric("CNV.focal")
})

#' @rdname CNV.focal
setMethod("CNV.focal", signature(object = "CNV.analysis"), function(object, sig_cgenes = FALSE, conf = 0.99, R = 100, blockLength = 500000, proportionLength = TRUE, ...){
  if(ncol(object@anno@genome) == 2) {
    stop("CNV.focal is not compatible with mouse arrays.")
  }

  if(!length(object@anno@detail) >= 1){
    stop("Please run CNV.detail() to specify genes of interest.")
  }

  amp_bins <- vector(mode='list', length=ncol(object@fit$ratio))
  del_bins <- vector(mode='list', length=ncol(object@fit$ratio))
  detail.regions.amp <- vector(mode='list', length=ncol(object@fit$ratio))
  detail.regions.del <- vector(mode='list', length=ncol(object@fit$ratio))
  cancer.genes.amp <- vector(mode='list', length=ncol(object@fit$ratio))
  cancer.genes.del <- vector(mode='list', length=ncol(object@fit$ratio))

  for(i in 1:ncol(object@fit$ratio)){
    message(paste(colnames(object@fit$ratio)[i]), " (",round(i/ncol(object@fit$ratio)*100, digits = 3), "%", ")", sep = "")

    # Perform k-means on segments to identify likely copy-number state
    segs <- object@seg$summary[[i]]
    segs <- GRanges(seqnames = segs$chrom, IRanges(start = segs$loc.start, end = segs$loc.end), seqinfo = Seqinfo(genome = "hg19"))
    seqlevels(segs) <- object@anno@genome$chr
    segs$seg.median <- object@seg$summary[[i]]$seg.median - object@bin$shift[i]
    segs$num.markers <- object@seg$summary[[i]]$num.mark
    segs <- sort(segs)

    segs.rep <- rep(segs$seg.median, segs$num.markers)
    km <- kmeans(segs.rep, centers=3, nstart=25)

    bins.log2 <- object@bin$ratio[[i]] - object@bin$shift[i]
    bins <- object@anno@bins[names(bins.log2)]
    bins$weight <- 1/object@bin$variance[[i]][names(bins.log2)]
    bins$log2 <- as.numeric(bins.log2)
    bins$state <- rank(km$centers[, 1])[km$cluster]
    seqinfo(bins) <- Seqinfo(genome = "hg19")

    seg <- lapply(seq_len(3), function(s){
      x <- reduce(bins[bins$state == s])
      mcols(x)$state <- s
      x
    })
    seg <- c(seg[[1]], seg[[2]], seg[[3]])

    # Bootstrap to identify thresholds for deletions and amplifications, for each copy-number state
    boots <- bootRanges(bins, blockLength = blockLength, R = R, seg = seg, exclude = object@anno@exclude, proportionLength = proportionLength)

    boots.state <- split(boots$log2, boots$state)
    boots.state.mean <- sapply(boots.state, mean)
    boots.state.sd <- sapply(boots.state, sd)
    boots.state.low <- mapply(qnorm, boots.state.mean, boots.state.sd, p=(1-conf)/2)    # assume normal distr, two-sided
    boots.state.high <- mapply(qnorm, boots.state.mean, boots.state.sd, p=1-(1-conf)/2)

    bins.total <- object@bin$ratio[[i]] - object@bin$shift[i]
    bins.amp <- sort(bins.total[bins.total >= boots.state.high[bins$state]], decreasing = TRUE)
    bins.del <- sort(bins.total[bins.total <= boots.state.low[bins$state]])

    amp_bins[[i]] <- bins.amp
    del_bins[[i]] <- bins.del

    # Detail regions are subjected to the dynamic thresholds
    detail.genes <- object@anno@detail
    d1 <- as.matrix(findOverlaps(query = detail.genes, subject = object@anno@probes))
    d1seg <- as.matrix(findOverlaps(query = detail.genes, subject = seg))
    d2 <- data.frame(gene = values(detail.genes)$name[d1[,"queryHits"]], probe = names(object@anno@probes[d1[, "subjectHits"]]), stringsAsFactors = FALSE)
    detail.genes.ratio <- sapply(split(object@fit$ratio[d2[, "probe"], i], d2[, "gene"]), median, na.rm = TRUE)[values(detail.genes)$name]
    detail.genes.ratio <- detail.genes.ratio - object@bin$shift[i]
    detail.genes.seg <- round(sapply(split(seg$state[d1seg[, "subjectHits"]], values(detail.genes)$name[d1seg[,"queryHits"]]), mean))[values(detail.genes)$name]

    detail.regions.amp[[i]] <- sort(detail.genes.ratio[which(detail.genes.ratio >= boots.state.high[detail.genes.seg])], decreasing = TRUE)
    detail.regions.del[[i]] <- sort(detail.genes.ratio[which(detail.genes.ratio <= boots.state.low[detail.genes.seg])])

    # Genes from Cancer Gene Census are subjected to the dynamic thresholds
    if(sig_cgenes) {

      if(object@anno@args$genome == "hg38"){
        data("consensus_cancer_genes_hg38")
      } else {
        data("consensus_cancer_genes_hg19")
      }

      cancer.genes <- cancer_genes[setdiff(cancer_genes$SYMBOL, object@anno@detail$name)]
      d1 <- as.matrix(findOverlaps(query = cancer.genes, subject = object@anno@probes))
      d1seg <- as.matrix(findOverlaps(query = cancer.genes, subject = seg))
      d2 <- data.frame(gene = values(cancer.genes)$SYMBOL[d1[, "queryHits"]], probe = names(object@anno@probes[d1[, "subjectHits"]]),stringsAsFactors = FALSE)
      cancer.genes.ratio <- sapply(split(object@fit$ratio[d2[, "probe"], i], d2[, "gene"]), median, na.rm = TRUE)[values(cancer.genes)$SYMBOL]
      cancer.genes.ratio <- cancer.genes.ratio - object@bin$shift[i]
      cancer.genes.seg <- round(sapply(split(seg$state[d1seg[, "subjectHits"]], values(cancer.genes)$SYMBOL[d1seg[,"queryHits"]]), mean, na.rm=TRUE))[values(cancer.genes)$SYMBOL]

      cancer.genes.amp[[i]] <- sort(cancer.genes.ratio[which(cancer.genes.ratio >= boots.state.high[cancer.genes.seg])], decreasing = TRUE)
      cancer.genes.del[[i]] <- sort(cancer.genes.ratio[which(cancer.genes.ratio <= boots.state.low[cancer.genes.seg])])
    }
  }

  names(amp_bins) <- colnames(object@fit$ratio)
  names(del_bins) <- colnames(object@fit$ratio)
  names(detail.regions.amp) <- colnames(object@fit$ratio)
  names(detail.regions.del) <- colnames(object@fit$ratio)
  names(cancer.genes.amp) <- colnames(object@fit$ratio)
  names(cancer.genes.del) <- colnames(object@fit$ratio)

  object@detail$amp.bins <- amp_bins
  object@detail$del.bins <- del_bins
  object@detail$amp.detail.regions <- detail.regions.amp
  object@detail$del.detail.regions <- detail.regions.del
  object@detail$amp.cancer.genes <- cancer.genes.amp
  object@detail$del.cancer.genes <- cancer.genes.del

  return(object)
})
