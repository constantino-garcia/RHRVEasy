#' @export
RHRVEasy <-
  function(folders,
           correction = TRUE,
           correctionMethod = "bonferroni",
           verbose = FALSE,
           format = "RR",
           typeAnalysis = c('fourier', 'wavelet'),
           significance_level = 0.05,
           nonLinear = FALSE,
           doRQA = FALSE,   # ignored if nonLinear = FALSE
           saveHRVindexesInPath = NULL,
           nJobs = 1,
           ...) {

    typeAnalysis <- match.arg(typeAnalysis)
    if (!correction) {
      correctionMethod = NULL
    }
    easyOptions <- buildEasyOptions(
      verbose = verbose,
      significance = significance_level,
      method = correctionMethod
    )
    HRVIndices <- calculateHRVIndices(
      folders = folders,
      format = format,
      typeAnalysis = typeAnalysis,
      nonLinear = nonLinear,
      doRQA = doRQA,
      nJobs = nJobs,
      easyOptions = easyOptions,
      ...
    )
    tempfilename <- tempfile(pattern = "RHRVEasy", tmpdir = tempdir(), fileext = ".RDS")
    saveRDS(HRVIndices, file = tempfilename)
    message(
      paste0("Saving temporary HRV Indices. You may read them by using readRDS('", tempfilename, "')")
    )

    pVals <- statsAnalysis(HRVIndices, easyOptions)
    results <- list("HRVIndices" = HRVIndices, "stats" = pVals)
    class(results) <- "RHRVEasyResult"
    attr(results, "easyOptions") <- easyOptions
    results
  }


#' @export
print.RHRVEasyResult <- function(x, digits = getOption("digits"), ...) {
  easyOptions <- attr(x, "easyOptions")
  significantx <- x$stats[x$stats$adj.p.value < easyOptions$significance, ]
  if (is.null(significantx) || nrow(significantx) == 0) {
    cat(
      "No significant differences were found between groups\n"
    )
  } else {
    junk <- foreach(test = iterators::iter(significantx, by = "row")) %do% {
      # Compute mean CIS using Normal CIs or Boostrapped CIs depending on the
      # data distribution
      grouped_data <- split(x$HRVIndices[[test$HRVIndex]], x$HRVIndices$group)
      if (test$method == "ANOVA") {
        cis <- lapply(grouped_data, \(x) {
          ci <- t.test(x)$conf.int
          attr(ci, "method") <- "Normal CI without adjustment"
          ci
        })
      } else {
        cis <- lapply(grouped_data, \(x) {
          conf <- 0.95
          ci <- boot::boot.ci(
            boot::boot(x, statistic = \(x, i) mean(x[i]), R = 1000),
            type = "basic",
            conf = conf
          )
          ci <- ci$basic[, 4:5]
          attr(ci, "conf.level") <- conf
          attr(ci, "method") <- "Boostrapped CI without adjustment"
          ci
        })
      }
      # Actual printing happens here
      with(test, {
        method <- ifelse(easyOptions$method == "none", "", easyOptions$method)
        cat(
          sep = "",
          "Significant differences in ",  HRVIndex, " (", method," p-value = ", format(adj.p.value, digits), "):\n"
        )
      })
      for (group in names(cis)) {
        conf <- attr(cis[[group]], "conf.level")
        method <- attr(cis[[group]], "method")
        cat(
          sep = "",
          "\t", group, "'s ", conf * 100, "% CI of the mean: (",
          round(cis[[group]][1], digits = digits), ", ",
          round(cis[[group]][2], digits = digits), ")",
          " [", method,"]\n"
        )
      }
      cat("\n")
    }
  }
  invisible(x)
}


buildEasyOptions <- function(verbose, significance, method) {
  # fake call to p.adjust to check for the method name
  invisible(p.adjust(rep(0.1, 3), method = method))
  stopifnot((significance > 0) && (significance < 1))
  list(
    "verbose" = verbose,
    "significance" = significance,
    "method" = method
  )
}


calculateHRVIndices <- function(
    folders,
    format,
    typeAnalysis,
    nonLinear,
    doRQA,   # ignored if nonLinear = FALSE
    nJobs,
    easyOptions,
    ...) {


    cl <- prepareEasyCluster(nJobs, easyOptions$verbose)
    easyOptions$parallel <- !is.null(cl)

    filesByFolder <- lapply(folders, \(folder) {
      fileValidation(folder, easyOptions)
      data.frame(
        "folder" = folder,
        "file" = list.files(folder),
        "group" = splitPath(folder)[[1]]
      )
    })
    files <- do.call(rbind, filesByFolder)
    files$group <- factor(files$group)

    timeResults <- easyTimeAnalysis(
      format = format,
      files = files$file,
      groups = files$group,
      paths = files$folder,
      easyOptions = easyOptions,
      ...
    )
    # Frequency analysis
    if (typeAnalysis == "fourier") {
      freqResults <- easyFreqAnalysis(
        format = format,
        files = files$file,
        groups = files$group,
        paths = files$folder,
        easyOptions = easyOptions,
        ...
      )
    }
    if (typeAnalysis == "wavelet") {
      freqResults <-
        easyWaveletAnalysis(
          format = format,
          files = files$file,
          groups = files$group,
          paths = files$folder,
          type = typeAnalysis,
          easyOptions = easyOptions,
          ...
        )
    }
    if (nonLinear) {
      nonlinearResults <-
        easyNonLinearAnalysis(
          format = format,
          files = files$file,
          groups = files$group,
          paths = files$folder,
          easyOptions = easyOptions,
          doRQA = doRQA,
          ...
        )
    }
    if (!is.null(cl)) {
      parallel::stopCluster(cl)
    }

    # Merge
    allResults <- merge(
      timeResults,
      freqResults,
      by = c("file", "group"),
      all = TRUE
    )
    if (nonLinear) {
      allResults <- merge(
        allResults,
        nonlinearResults,
        by = c("file", "group"),
        all = TRUE
      )
    }
    allResults
  }

