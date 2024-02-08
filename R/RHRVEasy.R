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
           nJobs = 1,
           saveHRVIndicesInPath = NULL,
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
    results <- RHRVEasyResult(HRVIndices, pVals, easyOptions)

    if (!is.null(saveHRVIndicesInPath)) {
      saveHRVIndices(results, saveHRVIndicesInPath)
    }

    results
  }

RHRVEasyResult <- function(HRVIndices, pVals, easyOptions) {
  results <- list("HRVIndices" = HRVIndices, "stats" = pVals)
  class(results) <- "RHRVEasyResult"
  attr(results, "easyOptions") <- easyOptions
  results
}

#' @export
RHRVEasyStats <- function(RHRVEasyResultObject,
                          correctionMethod = c("bonferroni", "holm", "hochberg", "hommel", "BH",
                                               "BY", "fdr", "none"),
                          significance = 0.05) {
  if (!inherits(RHRVEasyResultObject, "RHRVEasyResult")) {
    stop("RHRVEasyResultObject should be a 'RHRVEasyResult' object, as returned by 'RHRVEasy()'")
  }
  correctionMethod <- match.arg(correctionMethod)
  stopifnot((significance > 0) && (significance < 1))
  easyOptions <- attr(RHRVEasyResultObject, "easyOptions")
  easyOptions$method <- correctionMethod
  easyOptions$significance <- significance
  HRVIndices <- RHRVEasyResultObject$HRVIndices
  pVals <- statsAnalysis(HRVIndices, easyOptions)
  RHRVEasyResult(HRVIndices, pVals, easyOptions)
}



#' @importFrom writexl write_xlsx
#' @export
saveHRVIndices <- function(results, saveHRVIndicesInPath = ".") {
  tryCatch({
    filename <- file.path(
      saveHRVIndicesInPath,
      paste0(
        paste(unique(results$HRVIndices$group), collapse = "_Vs_"),
        ".xlsx"
      )
    )
    write_xlsx(results$HRVIndices, filename)
  },
  error = function(e) {
    message(paste0("Could not save indices in path '", saveHRVIndicesInPath, "'"))
    message(e)
  })
}



#' @importFrom boot boot.ci
#' @importFrom boot boot
computeEasyCIs <- function(easyObject, test, confLevel) {
  groupedData <- split(easyObject$HRVIndices[[test$HRVIndex]], easyObject$HRVIndices$group)
  if (test$method == "ANOVA") {
    cis <- lapply(groupedData, \(x) {
      ci <- t.test(x, conf.level = confLevel)$conf.int
      attr(ci, "method") <- "Normal CI without adjustment"
      ci
    })
  } else {
    cis <- lapply(groupedData, \(x) {
      ci <- boot::boot.ci(
        boot::boot(x[!is.na(x)], statistic = \(x, i) mean(x[i]), R = 1000),
        type = "basic",
        conf = confLevel
      )
      ci <- ci$basic[, 4:5]
      attr(ci, "conf.level") <- confLevel
      attr(ci, "method") <- "Bootstrap CI without adjustment"
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
    rep(" ", nspaces), group, "'s mean", conf * 100, "% CI: (",
    round(cis[[group]][1], digits = digits), ", ",
    round(cis[[group]][2], digits = digits), ")",
    " [", method,"]\n"
  )
}


#' @importFrom tibble as_tibble
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
        method <- test$pairwise[[1]]$method[[1]]
        isSignificantPosthoc <- test$pairwise[[1]]$adj.p.value < easyOptions$significance
        if (!any(isSignificantPosthoc)) {
          cat(
            rep(" ", firstLevelSpaces),
            "No significant differences were found between groups in post-hoc tests (",
            method, " + ", adjMethod, "-p-value adjustment).\n"
          )
        } else {
          cat(
            sep = "",
            rep(" ", firstLevelSpaces),
            "Significant differences in the post-hoc tests (",
            method, " + ", adjMethod, "-p-value adjustment):\n"
          )

          significantPosthocs <- test$pairwise[[1]][isSignificantPosthoc, ]
          significantPosthocs <- significantPosthocs[order(significantPosthocs$group1), ]
          posthocTable <- capture.output(
            as_tibble(significantPosthocs[, c("group1", "group2", "adj.p.value")])
          )
          # eliminate tibble head and type info
          posthocTable <- posthocTable[-c(1, 3)]
          maxLen <- max(sapply(posthocTable, nchar))
          printSpaces <- paste0(rep(" ", nPosthocSpaces), collapse="")
          posthocTable <- sapply(posthocTable, \(x) paste0(printSpaces, x))
          cat(paste0(posthocTable, collapse = "\n"))
          cat("\n", printSpaces, rep("-", maxLen), "\n", sep="")
          for (group in sort(names(cis))) {
            printGroupCI(cis, group, digits, nspaces = nPosthocSpaces)
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


