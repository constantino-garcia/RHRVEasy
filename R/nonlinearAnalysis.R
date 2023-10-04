# Select the range for the computation of nonlinear statistics
selectMaxEmbeddingDim <- function(embeddingDim) {embeddingDim + 2}

tryRQA <- function(hrv.data, file, config) {
  # TODO: easy.config
  message(paste("Doing RQA for ", file))

  hrv.data = tryCatch({
    hrv.data$NonLinearAnalysis[[1]]$rqa = with(
      config,
      nonlinearTseries::rqa(
        takens = NULL,
        time.series = hrv.data$Beat$RR,
        embedding.dim = embeddingDim,
        time.lag = timeLag,
        radius = smallCorrRadius,
        do.plot = FALSE,
        save.RM = FALSE
      )
    )
    gc()
    hrv.data
  },
  error = function(e) {
    warning(paste(
      "RQA failed for file",
      file
    ))
    hrv.data
  })
  hrv.data
}


findSmallRadius <- function(hrv.data, file, config) {
  cd = hrv.data$NonLinearAnalysis[[1]]$correlation$computations
  if (!is.null(cd)) {
    # An alternative strategy
    # correlations = cd$corr.matrix[cd$embedding.dims == config$embeddingDim,]
    correlations = colMeans(cd$corr.matrix)
    interpolatedFun = approxfun(cd$radius, correlations - config$kRefCorrelation)
    smallCorrRadius = tryCatch(
      uniroot(
        interpolatedFun, interval = range(cd$radius)
      )$root,
      error = function(e) {
        warning(
          paste(
            "Could not find a small radius Lyapunov/RQA in file",
            file,
            ". Using min(radius)"
          )
        )
        min(cd$radius)
      }
    )
  } else {
    warning(
      paste(
        "Could not find a small radius Lyanpunov/RQA in file",
        file,
        ". Using radius = 10"
      )
    )
    smallCorrRadius = 10 # TODO get from config
  }
  smallCorrRadius
}


collectRQA <- function(rqa) {
  RQAStats = c("REC", "RATIO", "DET", "DIV", "Lmax", "Lmean", "LmeanWithoutMain",
               "ENTR", "TREND", "LAM", "Vmax", "Vmean")
  # Set names so that lapply sets them in results
  names(RQAStats) = RQAStats
  resultsRQA = lapply(
    RQAStats,
    function(statName) {
      result = rqa[[statName]]
      isInvalid = function(result) {
        is.null(result) || is.infinite(result) || is.nan(result)
      }
      ifelse(isInvalid(result), NA, result)
    }
  )
  resultsRQA
}


collectResuls <- function(hrv.data, file, class, config, doRQA) {
  results = lapply(
    # apply an extraction function to these pairs of "Final name" = "nonlinear list name"
    c(
      "CorrelationStatistic" = "correlation",
      "SampleEntropy" = "sampleEntropy",
      "MaxLyapunov" = "lyapunov"
    ),
    function(analysisName) {
      nla = hrv.data$NonLinearAnalysis[[1]][[analysisName]]$statistic
      if (any(is.null(nla)) || all(is.na(nla))) {
        NA
      } else {
        mean(nla, na.rm = TRUE)
      }
    }
  )
  results$PoincareSD1 = hrv.data$NonLinearAnalysis[[1]]$PoincarePlot$SD1
  results$PoincareSD1 = ifelse(is.null(results$PoincareSD1), NA, results$PoincareSD1)
  results$PoincareSD2 = hrv.data$NonLinearAnalysis[[1]]$PoincarePlot$SD2
  results$PoincareSD2 = ifelse(is.null(results$PoincareSD2), NA, results$PoincareSD2)

  if (doRQA) {
    RQAResults = collectRQA(hrv.data$NonLinearAnalysis[[1]]$rqa)
    results = c(results, RQAResults)
  }

  configList = with(config,
                    list(
                      "embeddingDim" = embeddingDim,
                      "timeLag" = timeLag
                    )
  )
  results = c(
    list(
      "filename" = file,
      "group" = class
    ),
    results,
    configList
  )
  as.data.frame(results)
}

#' @importFrom RHRV CalculateTimeLag
tryTimeLagEstimation <- function(hrv.data, easy_options) {
  lag = 30
  timeLag = tryCatch({
    timeLag <-
      CalculateTimeLag(
        hrv.data,
        technique = "acf",
        method = "first.minimum",
        lagMax = lag,
        doPlot = FALSE
      )
    timeLag
  },
  error = function(cond) {
    tryCatch({
      timeLag <-
        CalculateTimeLag(
          hrv.data,
          technique = "acf",
          method = "first.e.decay",
          lagMax = lag,
          doPlot = FALSE
        )
      timeLag
    },
    error = function(cond) {
      tryCatch({
        timeLag <-
          CalculateTimeLag(
            hrv.data,
            technique = "ami",
            method = "first.minimum",
            lagMax = lag,
            doPlot = FALSE
          )
        timeLag
      },
      error = function(cond) {
        tryCatch({
          timeLag <-
            CalculateTimeLag(
              hrv.data,
              technique = "ami",
              method = "first.e.decay",
              lagMax = lag,
              doPlot = FALSE
            )
          timeLag
        },
        error = function(cond) {
          if (easy_options$verbose) {
            message("Using default timeLag for current recording...")
          }
          30
        })
      })
    })
  })
  if (easy_options$verbose) {
    message(c("Time Lag for takens reconstruction: ", timeLag))
  }
  timeLag
}

tryEmbeddingDim <- function(hrv.data, file, config, embeddingDimDefault = 15) {
  embeddingDim = tryCatch(
    with(config,
         CalculateEmbeddingDim(
           hrv.data,
           numberPoints = 10000,
           timeLag = timeLag,
           maxEmbeddingDim = 15,
           threshold = 0.90,
           doPlot = showNonLinearPlots
         )
    ),
    error = function(e) {
      NA
    }
  )
  # TODO: unifiy is.na with 0
  if (is.na(embeddingDim) || (embeddingDim == 0)) {
    embeddingDim = embeddingDimDefault
    warning(paste(
      "Proper embedding dim not found for file",
      file,
      "Setting to 15."
    ))
  }
  embeddingDim
}

augmentCorrDimCalculations <- function(hrv.data, increaseRadius, config) {
  hrv.data = with(config, {
    if (increaseRadius) {
      minRadius = kMaxRadius + 10
      maxRadius = max(kMaxRadius + 20, 2 * kMaxRadius)
      pointsRadius = 10
    } else {
      minRadius = 1
      maxRadius = 9
      pointsRadius = 5
    }
    hrv.data = CreateNonLinearAnalysis(hrv.data)
    CalculateCorrDim(
      hrv.data,
      indexNonLinearAnalysis = 2,
      minEmbeddingDim = embeddingDim,
      maxEmbeddingDim = maxEmbeddingDim,
      timeLag = timeLag,
      minRadius = minRadius,
      maxRadius = maxRadius,
      pointsRadius = pointsRadius,
      theilerWindow = kTheilerWindow,
      corrOrder = 2,
      doPlot = showNonLinearPlots
    )
  })
  cd1Computations = hrv.data$NonLinearAnalysis[[1]]$correlation$computations
  cd2Computations = hrv.data$NonLinearAnalysis[[2]]$correlation$computations
  if (increaseRadius) {
    cd1Computations$corr.matrix = cbind(cd2Computations$corr.matrix,
                                        cd1Computations$corr.matrix)
  } else {
    cd1Computations$corr.matrix = cbind(cd1Computations$corr.matrix,
                                        cd2Computations$corr.matrix)
  }
  cd1Computations$radius = as.numeric(colnames(cd1Computations$corr.matrix))
  # radius should be monotonic decreasing
  stopifnot(all(diff(cd1Computations$radius) < 0))
  hrv.data$NonLinearAnalysis[[1]]$correlation$computations = cd1Computations
  hrv.data$NonLinearAnalysis[[2]] = NULL
  hrv.data
}

tryCorrDimCalculations <- function(hrv.data, file, config) {
  hrv.data = tryCatch({
    hrv.data = with(
      config,
      CalculateCorrDim(
        hrv.data,
        indexNonLinearAnalysis = 1,
        minEmbeddingDim = embeddingDim,
        maxEmbeddingDim = maxEmbeddingDim,
        timeLag = timeLag,
        minRadius = config$kMinRadius,
        maxRadius = config$kMaxRadius,
        pointsRadius = 20,
        theilerWindow = kTheilerWindow,
        corrOrder = 2,
        doPlot = showNonLinearPlots
      )
    )
    cd1Computations = hrv.data$NonLinearAnalysis[[1]]$correlation$computations
    meanCorrMatrix = colMeans(cd1Computations$corr.matrix)
    # Increase radius to find small correlation values if needed
    if (max(meanCorrMatrix) < 0.25) {
      hrv.data = augmentCorrDimCalculations(hrv.data, increaseRadius = TRUE, config = config)
    }
    # Decrease radius to find small correlation values if needed
    if ((min(meanCorrMatrix) > config$kRefCorrelation) &&
        (min(cd1Computations$radius) == config$kMinRadius)) {
      hrv.data = augmentCorrDimCalculations(hrv.data, increaseRadius = FALSE, config = config)
    }
    hrv.data
  },
  error = function(error) {
    warning(
      paste(
        "Correlation Dimension calculations failed for file",
        file
      )
    )
    hrv.data
  })
  hrv.data
}

tryCorrDimEstimation <- function(hrv.data, file, config) {
  hrv.data = tryCatch({
    cd = hrv.data$NonLinearAnalysis[[1]]$correlation$computations
    stopifnot(!is.null(cd))

    filteredCd = nltsFilter(cd, threshold = 0.99)
    cdScalingRegion =
      estimate_scaling_region(filteredCd,
                              numberOfLinearRegions = 3,
                              doPlot = config$showNonLinearPlots)
    if (!cdScalingRegion$reliable) {
      warning(
        paste(
          "Scaling Region for file",
          file,
          "is not reliable.",
          "CorrDim and  SampleEntropy statistics may be wrong"
        )
      )
    }
    cdScalingRegion = cdScalingRegion$scalingRegion
    # nbPoints in scaling region should be at least two
    nbPoints = sum((cd$radius >= cdScalingRegion[1]) & (cd$radius <= cdScalingRegion[2]))
    if (nbPoints < 2) {
      # increase scalingRegion to ensure two points
      low = cd$radius[cd$radius < cdScalingRegion[1]]
      high = cd$radius[cd$radius > cdScalingRegion[2]]
      if (length(low) > 0) {
        cdScalingRegion[1] = max(low)
      }
      if (length(high) > 0) {
        cdScalingRegion[2] = min(high)
      }
      nbPoints = sum((cd$radius >= cdScalingRegion[1]) & (cd$radius <= cdScalingRegion[2]))
      # if the previous fails, increase the size of the scaling region slowly until meeting the criteria
      while (nbPoints < 2) {
        cdScalingRegion = cdScalingRegion * c(0.95, 1.05)
        nbPoints = sum((cd$radius >= cdScalingRegion[1]) & (cd$radius <= cdScalingRegion[2]))
      }
    }
    hrv.data = with(
      config,
      EstimateCorrDim(
        hrv.data,
        indexNonLinearAnalysis = 1,
        regressionRange = cdScalingRegion,
        useEmbeddings = embeddingDim:maxEmbeddingDim,
        doPlot = showNonLinearPlots
      )
    )
    attr(hrv.data$NonLinearAnalysis[[1]]$correlation$statistic, "automaticScalingRegion") = cdScalingRegion
    hrv.data
  },
  error = function(e) {
    warning(
      paste(
        "Correlation Dimension estimation failed for file",
        file,
        "Setting to NA"
      )
    )
    hrv.data$NonLinearAnalysis[[1]]$correlation$statistic = NA # TODO
    hrv.data
  })
  hrv.data
}

tryCorrDim <- function(hrv.data, file, config) {
  hrv.data = tryCorrDimCalculations(hrv.data, file, config)
  hrv.data = tryCorrDimEstimation(hrv.data, file, config)
  hrv.data
}


trySampleEntropy <- function(hrv.data, file, config) {
  hrv.data = trySampleEntropyCalculations(hrv.data, file, config)
  hrv.data = trySampleEntropyEstimation(hrv.data, file, config)
  hrv.data
}

trySampleEntropyCalculations <- function(hrv.data, file, config) {
  hrv.data = tryCatch({
    stopifnot(!is.null(hrv.data$NonLinearAnalysis[[1]]$correlation$computations))
    with(
      config,
      CalculateSampleEntropy(hrv.data,
                             indexNonLinearAnalysis = 1,
                             doPlot = config$showNonLinearPlots)
    )
  },
  error = function(error) {
    warning(
      paste(
        "Sample Entropy calculations failed for file",
        file
      )
    )
    hrv.data
  })
  hrv.data
}


# TODO: use the attribute set in correlation dimension
trySampleEntropyEstimation <- function(hrv.data, file, config) {
  hrv.data = tryCatch({
    EstimateSampleEntropy(hrv.data,
                          indexNonLinearAnalysis = 1,
                          doPlot = config$showNonLinearPlots)
  },
  error = function(e) {
    warning(
      paste(
        "Sample entropy estimation failed for file",
        file,
        "Setting to NA"
      )
    )
    hrv.data$NonLinearAnalysis[[1]]$sampleEntropy$statistic = NA
    hrv.data
  })
  hrv.data
}

tryLyapunov <- function(hrv.data, file, config) {
  hrv.data = tryLyapunovCalculations(hrv.data, file, config)
  hrv.data = tryLyapunovEstimation(hrv.data, file, config)
  hrv.data
}

tryLyapunovCalculations <- function(hrv.data, file, config) {
    hrv.data = tryCatch(
      with(config,
           CalculateMaxLyapunov(
             hrv.data,
             indexNonLinearAnalysis = 1,
             minEmbeddingDim = embeddingDim,
             maxEmbeddingDim = maxEmbeddingDim,
             timeLag = timeLag,
             radius = smallCorrRadius,
             theilerWindow = kTheilerWindow,
             doPlot = showNonLinearPlots
           )
      ),
      error = function(e) {
        warning(
          paste(
            "Lyapunov computations failed for file",
            file
          )
        )
        hrv.data
      }
    )
    hrv.data
}

tryLyapunovEstimation <- function(hrv.data, file, config) {
  hrv.data = tryCatch({
    lyapComp = hrv.data$NonLinearAnalysis[[1]]$lyapunov$computations
    stopifnot(!is.null(lyapComp))
    lyapunovScalingRegion = estimate_scaling_region(lyapComp)
    hrv.data = with(
      config,
      EstimateMaxLyapunov(
        hrv.data,
        indexNonLinearAnalysis = 1,
        regressionRange = lyapunovScalingRegion,
        useEmbeddings = embeddingDim:maxEmbeddingDim,
        doPlot = showNonLinearPlots
      )
    )
  },
  error = function(error) {
    warning(
      paste0(
        "Lyapunov Estimation failed for file ",
        file,
        ". Setting to NA."
      )
    )
    hrv.data$NonLinearAnalysis[[1]]$lyapunov$statistic = NA
    hrv.data
  })
  hrv.data
}


nlaSingleFile <- function(file, rrs2, format, class, easy_options, doRQA) {
  # Config is a list that sets some parameters needed for the analysis
  config = list(
    # Set the value of the Theiler window for avoiding temporal correlations
    "kTheilerWindow" = 10,
    # Correlation value used to get a reasonable radious for both RQA and Lyapunov
    "kRefCorrelation" = 1e-3,
    # Min and max radius for correlation dimension
    "kMinRadius" = 10,
    "kMaxRadius" = 50
  )

  hrv.data = preparing_analysis(file = file, rrs = rrs2, format = format,
                                easy_options = easy_options)
  hrv.data = CreateNonLinearAnalysis(hrv.data)
  config$timeLag = tryTimeLagEstimation(hrv.data, easy_options)

  #Poincare does not depend on the calculation of time lag or correlation dimension
  #unlike the rest of the nonlinear statistics, its calculation should never fail
  hrv.data = PoincarePlot(hrv.data, indexNonLinearAnalysis = 1, timeLag = 1)

  #Set to TRUE to display correlation dimension calculation and Lyapunov related plots
  config$showNonLinearPlots = FALSE

  if (config$showNonLinearPlots) {
    PlotNIHR(hrv.data, main = paste("NIHR of ", file))
  }

  config$embeddingDim = tryEmbeddingDim(hrv.data, file, config, embeddingDimDefault = 15)
  config$maxEmbeddingDim = selectMaxEmbeddingDim(config$embeddingDim)

  hrv.data = NonLinearNoiseReduction(HRVData = hrv.data,
                                     embeddingDim = config$embeddingDim,
                                     radius = NULL)


  hrv.data = tryCorrDim(hrv.data, file, config)
  hrv.data = trySampleEntropy(hrv.data, file, config)

  ## Get a reasonable radius for both Lyapunov and RQA
  config$smallCorrRadius = findSmallRadius(hrv.data, file, config)
  hrv.data = tryLyapunov(hrv.data, file, config)
  if (doRQA) {
    hrv.data = tryRQA(hrv.data, file, config)
  }

  results = collectResuls(hrv.data, file, class, config, doRQA)
  if (easy_options$verbose) {
    message(
      paste(
        capture.output({
          print(results)
        }),
        collapse = "\n"
      )
    )
  }
  results
}



#' @importFrom foreach foreach
#' @importFrom RHRV CreateNonLinearAnalysis CalculateEmbeddingDim
#' @importFrom RHRV PlotNIHR NonLinearNoiseReduction
#' @importFrom RHRV CalculateCorrDim EstimateCorrDim
#' @importFrom RHRV CalculateSampleEntropy EstimateSampleEntropy
#' @importFrom RHRV CalculateMaxLyapunov EstimateMaxLyapunov
#' @importFrom RHRV RQA PoincarePlot
#' @importFrom nonlinearTseries rqa
non_linear_analysis <-
  function(format, files, class, rrs2, easy_options, doRQA, ...) {
    # Compute all nonlinear statistics but RQA in parallel using %dopar%.
    # RQA is avoided due to its high memory consumption
    resultsDataFrame =
      foreach(
      file = files,
      .combine = rbind.data.frame,
      # .export = c("preparing_analysis", "easy_call",
      #             "nltsFilter", "selectMaxEmbeddingDim", "tryRQA",
      #             "findSmallRadius", "collectRQA", "collectResuls",
      #             "tryTimeLagEstimation", "tryEmbeddingDim", "augmentCorrDimCalculations",
      #             "tryCorrDimCalculations", "tryCorrDimEstimation", "tryCorrDim",
      #             "trySampleEntropy", "trySampleEntropyCalculations",
      #             "trySampleEntropyEstimation", "tryLyapunov",
      #             "tryLyapunovCalculations", "tryLyapunovEstimation",
      #             "nlaSingleFile"
      # ),
      # .packages = "RHRV",
      .errorhandling = "pass"
    ) %dopar% {
      suppressPackageStartupMessages(library("RHRV", character.only = TRUE))
      fileResults = nlaSingleFile(file, rrs2, format, class, easy_options, doRQA)
      fileResults
    } # end of %dopar%
    resultsDataFrame
  }


# Dunn Statistical tests for non-linear statistics
dunnNonLinear <- function(dfM, correctionMethod, easy_options = easy_options) {
  dfM$group = factor(dfM$group)
  CorrelationStatistic = NA
  SampleEntropy = NA
  MaxLyapunov  = NA
  REC = NA
  RATIO = NA
  DET = NA
  DIV = NA
  Lmax = NA
  Lmean = NA
  LmeanWithoutMain = NA
  ENTR = NA
  TREND = NA
  LAM = NA
  Vmax = NA
  Vmean = NA
  PoincareSD1 = NA
  PoincareSD2  = NA


  CorrelationStatistic = posthoc.kruskal.dunn.test.CheckAllValuesEqual(CorrelationStatistic ~ group,
                                                                       data = dfM,
                                                                       p.adjt = correctionMethod, easy_options = easy_options)
  SampleEntropy = posthoc.kruskal.dunn.test.CheckAllValuesEqual(SampleEntropy ~ group, data =
                                                                  dfM, p.adjt = correctionMethod, easy_options = easy_options)
  MaxLyapunov = posthoc.kruskal.dunn.test.CheckAllValuesEqual(MaxLyapunov ~ group, data =
                                                                dfM, p.adjt = correctionMethod, easy_options = easy_options)
  REC = posthoc.kruskal.dunn.test.CheckAllValuesEqual(REC ~ group, data =
                                                        dfM, p.adjt = correctionMethod, easy_options = easy_options)
  RATIO = posthoc.kruskal.dunn.test.CheckAllValuesEqual(RATIO ~ group, data =
                                                          dfM, p.adjt = correctionMethod, easy_options = easy_options)
  DET = posthoc.kruskal.dunn.test.CheckAllValuesEqual(DET ~ group, data =
                                                        dfM, p.adjt = correctionMethod, easy_options = easy_options)
  DIV = posthoc.kruskal.dunn.test.CheckAllValuesEqual(DIV ~ group, data =
                                                        dfM, p.adjt = correctionMethod, easy_options = easy_options)
  Lmax = posthoc.kruskal.dunn.test.CheckAllValuesEqual(Lmax ~ group, data =
                                                         dfM, p.adjt = correctionMethod, easy_options = easy_options)
  Lmean = posthoc.kruskal.dunn.test.CheckAllValuesEqual(Lmean ~ group, data =
                                                          dfM, p.adjt = correctionMethod, easy_options = easy_options)
  LmeanWithoutMain = posthoc.kruskal.dunn.test.CheckAllValuesEqual(LmeanWithoutMain ~ group, data =
                                                                     dfM, p.adjt = correctionMethod, easy_options = easy_options)
  ENTR = posthoc.kruskal.dunn.test.CheckAllValuesEqual(ENTR ~ group, data =
                                                         dfM, p.adjt = correctionMethod, easy_options = easy_options)
  TREND = posthoc.kruskal.dunn.test.CheckAllValuesEqual(TREND ~ group, data =
                                                          dfM, p.adjt = correctionMethod, easy_options = easy_options)
  LAM = posthoc.kruskal.dunn.test.CheckAllValuesEqual(LAM ~ group, data =
                                                        dfM, p.adjt = correctionMethod, easy_options = easy_options)
  Vmax = posthoc.kruskal.dunn.test.CheckAllValuesEqual(Vmax ~ group, data =
                                                         dfM, p.adjt = correctionMethod, easy_options = easy_options)
  Vmean = posthoc.kruskal.dunn.test.CheckAllValuesEqual(Vmean ~ group, data =
                                                          dfM, p.adjt = correctionMethod, easy_options = easy_options)
  PoincareSD1 = posthoc.kruskal.dunn.test.CheckAllValuesEqual(PoincareSD1 ~ group, data =
                                                                dfM, p.adjt = correctionMethod, easy_options = easy_options)
  PoincareSD2  = posthoc.kruskal.dunn.test.CheckAllValuesEqual(PoincareSD2 ~ group, data =
                                                                 dfM, p.adjt = correctionMethod, easy_options = easy_options)

  list(
    CorrelationStatistic,
    SampleEntropy,
    MaxLyapunov,
    REC,
    RATIO,
    DET,
    DIV,
    Lmax,
    Lmean,
    LmeanWithoutMain,
    ENTR,
    TREND,
    LAM,
    Vmax,
    Vmean,
    PoincareSD1,
    PoincareSD2
  )
}
