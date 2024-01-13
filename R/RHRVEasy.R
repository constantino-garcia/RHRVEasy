#' @export
RHRVEasy <-
  function(folders,
           correctionMethod = c("bonferroni", "holm", "hochberg", "hommel", "BH",
                                "BY", "fdr", "none"),
           verbose = FALSE,
           format = "RR",
           typeAnalysis = c('fourier', 'wavelet'),
           significance = 0.05,
           nonLinear = FALSE,
           doRQA = FALSE,   # ignored if nonLinear = FALSE
           saveHRVindexesInPath = NULL,
           nJobs = 1,
           ...) {

    typeAnalysis <- match.arg(typeAnalysis)
    correctionMethod <- match.arg(correctionMethod)
    easyOptions <- buildEasyOptions(
      verbose = verbose,
      significance = significance,
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


computeEasyCIs <- function(easyObject, test, confLevel) {
  grouped_data <- split(easyObject$HRVIndices[[test$HRVIndex]], easyObject$HRVIndices$group)
  if (test$method == "ANOVA") {
    cis <- lapply(grouped_data, \(x) {
      ci <- t.test(x, conf.level = confLevel)$conf.int
      attr(ci, "method") <- "Normal CI without adjustment"
      ci
    })
  } else {
    cis <- lapply(grouped_data, \(x) {
      ci <- boot::boot.ci(
        boot::boot(x[!is.na(x)], statistic = \(x, i) mean(x[i]), R = 1000),
        type = "basic",
        conf = confLevel
      )
      ci <- ci$basic[, 4:5]
      attr(ci, "conf.level") <- confLevel
      attr(ci, "method") <- "Boostrapped CI without adjustment"
      ci
    })
  }
  cis
}

printGroupCI <- function(cis, group, digits, nspaces = 2) {
  conf <- attr(cis[[group]], "conf.level")
  method <- attr(cis[[group]], "method")
  cat(
    sep = "",
    rep(" ", nspaces), group, "'s ", conf * 100, "% CI of the mean: (",
    round(cis[[group]][1], digits = digits), ", ",
    round(cis[[group]][2], digits = digits), ")",
    " [", method,"]\n"
  )
}

#' @export
print.RHRVEasyResult <- function(x, digits = getOption("digits"), ...) {
  firstLevelSpaces <- 2
  nPosthocSpaces <- 4
  easyOptions <- attr(x, "easyOptions")
  significantx <- x$stats[x$stats$adj.p.value < easyOptions$significance, ]
  if (is.null(significantx) || nrow(significantx) == 0) {
    cat(
      "No significant differences were found between groups\n"
    )
  } else {
    adjMethod <- ifelse(easyOptions$method == "none", "", easyOptions$method)
    for (sigRow in seq_len(nrow(significantx))) {
      test <- significantx[sigRow, ]
      # Compute mean CIS using Normal CIs or Boostrapped CIs depending on the
      # data distribution
      print(easyOptions$significance)
      cis <- computeEasyCIs(x, test, confLevel = 1 - easyOptions$significance)
      # Actual printing happens here
      with(test, {
        cat(
          sep = "",
          "Significant differences in ",  HRVIndex, " (", method, ", ",
          adjMethod," p-value = ", format(adj.p.value, digits), "):\n"
        )
      })
      if (length(cis) == 2) {
        for (group in names(cis)) {
          printGroupCI(cis, group, digits, nspaces = firstLevelSpaces)
        }
      } else {
        isSignificantPosthoc <- test$pairwise[[1]]$adj.p.value < easyOptions$significance
        if (!any(isSignificantPosthoc)) {
          cat(
            rep(" ", firstLevelSpaces),
            "No significant differences were found between groups in post-hoc tests\n"
          )
        } else {
          significantPosthocs <- test$pairwise[[1]][isSignificantPosthoc, ]
          for (prow in seq_along(significantPosthocs$group1)) {
            group1 <- significantPosthocs[prow, "group1"][[1]]
            group2 <- significantPosthocs[prow, "group2"][[1]]
            method <- significantPosthocs[prow, "method"][[1]]
            adj.p.value <- significantPosthocs[prow, "adj.p.value"][[1]]
            cat(
              sep = "",
              rep(" ", firstLevelSpaces),
              "Significant differences in the post-hoc comparison of ",
                group1, " and ", group2,
                " (", method, ", ", adjMethod, " p-value = ", format(adj.p.value, digits), "):\n"
            )
            printGroupCI(cis, group1, digits, nspaces = nPosthocSpaces)
            printGroupCI(cis, group2, digits, nspaces = nPosthocSpaces)
          }
        }
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

