# Creating time analysis data frames
#' @importFrom RHRV CreateTimeAnalysis
#' @importFrom foreach foreach
#' @importFrom foreach %dopar% %do%
easyTimeAnalysis <-
  function(format, files, groups, paths, easyOptions, ...) {
    opts <- NULL
    if (easyOptions$verbose) {
      opts <- list("progress" = updateProgressFactory("Time analysis", files))
    }
    dataFrame <- suppressPackageStartupMessages({
      foreach(
        file = files,
        itcounter = seq_along(files),
        group = groups,
        path = paths,
        .combine = rbind.data.frame,
        #.export = c("prepareAnalysis", "easyCall"),
        .errorhandling = "pass",
        .options.snow = opts
      ) %dopar% {
        suppressWarnings(
          suppressPackageStartupMessages(library("RHRV", character.only = TRUE, warn.conflicts = FALSE))
        )
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
        if (easyOptions$verbose && !easyOptions$parallel) {
          opts$progress(itcounter)
        }
        as.data.frame(row_list)
      }
    })
    dataFrame
  }