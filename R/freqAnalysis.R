
#' @importFrom foreach foreach
#' @importFrom RHRV CreateFreqAnalysis
#' @importFrom RHRV InterpolateNIHR CalculatePSD CalculateEnergyInPSDBands
easyFreqAnalysis <-
  function(format, files, groups, paths, easyOptions, ...) {
    opts <- NULL
    if (easyOptions$verbose) {
      opts <- list("progress" = updateProgressFactory("Freq analysis", files))
    }
    dataFrame <- foreach(
      file = files,
      itcounter = seq_along(files),
      group = groups,
      path = paths,
      .combine = rbind.data.frame,
      .export = c("prepareAnalysis", "easyCall"),
      # .packages = "RHRV",
      .errorhandling = "pass",
      .options.snow = opts
    ) %dopar% {
      suppressWarnings(
        suppressPackageStartupMessages(
          library("RHRV", character.only = TRUE, warn.conflicts = FALSE)
        )
      )
      hrv.data <- prepareAnalysis(file = file, rrs = path, format = format,
                                    easyOptions = easyOptions)
      hrv.data <- withCallingHandlers(
        {
          easyCall(hrv.data, InterpolateNIHR, ...)
        },
        warning = function(w) {
          w$message <- paste0("File ", file,": ", w$message, "\n")
          warning(w)
          invokeRestart("muffleWarning")
        }
      )
      zero_indexes <- which(hrv.data$HR == 0)
      hr_median <- median(hrv.data$HR[-zero_indexes])
      hrv.data$HR[zero_indexes] <- hr_median
      hrv.data <- easyCall(hrv.data, CreateFreqAnalysis, ...)
      hrv.data <- easyCall(hrv.data, CalculatePSD, doPlot = F, ...)
      x1 <- easyCall(hrv.data, CalculateEnergyInPSDBands, ...)
      names(x1) <- c("ULF", "VLF", "LF", "HF")
      row_list <- c(list("file" = file), x1, list("group" = group))
      if (easyOptions$verbose && !easyOptions$parallel) {
        opts$progress(itcounter)
      }
      as.data.frame(row_list)
    }
    dataFrame
  }

#' @importFrom foreach foreach
#' @importFrom RHRV CreateFreqAnalysis
#' @importFrom RHRV InterpolateNIHR CalculatePowerBand
easyWaveletAnalysis <-
  function(format, files, groups, paths, easyOptions, ...) {
    opts <- NULL
    if (easyOptions$verbose) {
      opts <- list("progress" = updateProgressFactory("Wavelet analysis", files))
    }
    freqResults <- foreach(
      file = files,
      itcounter = seq_along(files),
      group = groups,
      path = paths,
      .combine = rbind.data.frame,
      .export = c("prepareAnalysis", "easyCall"),
      .errorhandling = "pass",
      .options.snow = opts
    ) %dopar% {
      suppressWarnings(
        suppressPackageStartupMessages(
          library("RHRV", character.only = TRUE, warn.conflicts = FALSE)
        )
      )
      hrv.data <- prepareAnalysis(file = file, rrs = path, format = format, easyOptions = easyOptions)
      hrv.data <- withCallingHandlers(
        {
          easyCall(hrv.data, InterpolateNIHR, ...)
        },
        warning = function(w) {
          w$message <- paste0("File ", file,": ", w$message, "\n")
          warning(w)
          invokeRestart("muffleWarning")
        }
      )
      zero_indexes <- which(hrv.data$HR == 0)
      hr_median <- median(hrv.data$HR[-zero_indexes])
      hrv.data$HR[zero_indexes] <- hr_median

      hrv.data <- easyCall(hrv.data, CreateFreqAnalysis, ...)
      hrv.data <- SetVerbose(hrv.data, FALSE) # Set to False to avoid clutter
      hrv.data <- easyCall(hrv.data, CalculatePowerBand, ...)

      index <- length(hrv.data$FreqAnalysis)
      resultsWavelet <- hrv.data$FreqAnalysis[[index]]
      resultsWavelet$file <- file
      resultsWavelet$ULF <- sum(hrv.data$FreqAnalysis[[index]]$ULF)
      resultsWavelet$VLF <- sum(hrv.data$FreqAnalysis[[index]]$VLF)
      resultsWavelet$LF <- sum(hrv.data$FreqAnalysis[[index]]$LF)
      resultsWavelet$HF <- sum(hrv.data$FreqAnalysis[[index]]$HF)
      resultsWavelet[
        c("HRV", "LFHF", "Time", "wavelet", "bandtolerance", "depth", "type",
          paste0(c("ULF", "VLF", "LF", "HF"), rep(c("min", "max"), each = 4))
        )
      ] <- NULL
      row_list <- c(resultsWavelet, list("group" = group))
      if (easyOptions$verbose && !easyOptions$parallel) {
        opts$progress(itcounter)
      }
      as.data.frame(row_list)
    }
    freqResults
  }
