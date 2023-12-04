
#' @importFrom foreach foreach
#' @importFrom RHRV CreateFreqAnalysis
#' @importFrom RHRV InterpolateNIHR CalculatePSD CalculateEnergyInPSDBands
easyFreqAnalysis <-
  function(format, files, groups, paths, easyOptions, ...) {
    dataFrame <- foreach(
      file = files,
      group = groups,
      path = paths,
      .combine = rbind.data.frame,
      .export = c("prepareAnalysis", "easyCall"),
      # .packages = "RHRV",
      .errorhandling = "pass"
    ) %dopar% {
      suppressPackageStartupMessages(library("RHRV", character.only = TRUE, warn.conflicts = FALSE))
      hrv.data <- prepareAnalysis(file = file, rrs = path, format = format,
                                    easyOptions = easyOptions)
      hrv.data <- withCallingHandlers(
        {
          easyCall(hrv.data, InterpolateNIHR, ...)
        },
        warning = function(w) {
          w$message <- paste0("File ", file,": ", w$message)
          warning(w)
          invokeRestart("muffleWarning")
        }
      )
      zero_indexes <- which(hrv.data$HR == 0)
      hr_median <- median(hrv.data$HR[-zero_indexes])
      hrv.data$HR[zero_indexes] <- hr_median
      hrv.data <- easyCall(hrv.data, CreateFreqAnalysis, ...)
      hrv.data <- easyCall(hrv.data, CalculatePSD, doPlot = F, ...)
      # TODO: simplify
      x1 <- easyCall(hrv.data, CalculateEnergyInPSDBands, ...)
      names(x1) <- c("ULF", "VLF", "LF", "HF")
      row_list <- c(list("file" = file), x1, list("group" = group))
      as.data.frame(row_list)
    }
    dataFrame
  }

#' @importFrom foreach foreach
#' @importFrom RHRV CreateFreqAnalysis
#' @importFrom RHRV InterpolateNIHR CalculatePowerBand
easyWaveletAnalysis <-
  function(format, files, groups, paths, easyOptions, ...) {
    freqResults <- foreach(
      file = files,
      group = groups,
      path = paths,
      .combine = rbind.data.frame,
      .export = c("prepareAnalysis", "easyCall"),
      # .packages = "RHRV",
      .errorhandling = "pass"
    ) %dopar% {
      suppressPackageStartupMessages(library("RHRV", character.only = TRUE, warn.conflicts = FALSE))
      hrv.data <- prepareAnalysis(file = file, rrs = path, format = format, easyOptions = easyOptions)
      hrv.data <- withCallingHandlers(
        {
          easyCall(hrv.data, InterpolateNIHR, ...)
        },
        warning = function(w) {
          w$message <- paste0("File ", file,": ", w$message)
          warning(w)
          invokeRestart("muffleWarning")
        }
      )
      zero_indexes <- which(hrv.data$HR == 0)
      hr_median <- median(hrv.data$HR[-zero_indexes])
      hrv.data$HR[zero_indexes] <- hr_median

      hrv.data <- easyCall(hrv.data, CreateFreqAnalysis, ...)
      hrv.data <- SetVerbose(hrv.data, easyOptions$verbose)
      hrv.data <- easyCall(hrv.data, CalculatePowerBand, ...)

      index <- length(hrv.data$FreqAnalysis)
      resultsWavelet <- hrv.data$FreqAnalysis[[index]]
      resultsWavelet$file <- file
      resultsWavelet$HRV <- NULL
      resultsWavelet$ULF <- sum(hrv.data$FreqAnalysis[[index]]$ULF)
      resultsWavelet$VLF <- sum(hrv.data$FreqAnalysis[[index]]$VLF)
      resultsWavelet$LF <- sum(hrv.data$FreqAnalysis[[index]]$LF)
      resultsWavelet$HF <- sum(hrv.data$FreqAnalysis[[index]]$HF)
      resultsWavelet$LFHF <- NULL
      resultsWavelet$Time <- NULL
      row_list <- c(resultsWavelet, list("group" = group))
      as.data.frame(row_list)
    }
    freqResults
  }
