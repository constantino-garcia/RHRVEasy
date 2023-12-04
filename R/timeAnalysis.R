# Creating time analysis data frames
#' @importFrom RHRV CreateTimeAnalysis
#' @importFrom foreach foreach
#' @importFrom foreach %dopar% %do%
easyTimeAnalysis <-
  function(format, files, groups, paths, easyOptions, ...) {
    dataFrame <- foreach(
      file = files,
      group = groups,
      path = paths,
      .combine = rbind.data.frame,
      .errorhandling = "pass"
    ) %dopar% {
      suppressPackageStartupMessages(library("RHRV", character.only = TRUE, warn.conflicts = FALSE))
      hrv.data <- prepareAnalysis(file = file, rrs = path, format = format,
                                    easyOptions = easyOptions)
      hrv.data <- easyCall(hrv.data, CreateTimeAnalysis, ...)
      results <- hrv.data$TimeAnalysis[[1]]
      results$size <- NULL
      row_list <- c(
        list("file" = file),
        list("group" = group),
        results
      )
      as.data.frame(row_list)
    }
    dataFrame
  }
