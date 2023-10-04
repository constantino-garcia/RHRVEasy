file_validation <- function(path) {
  # 1. Check if path really exists
  if (dir.exists(path) != TRUE) {
    stop("\nThe path ", path, " does not exist")
  } else{
    cat("\nThe path ", path, " exists ")
  }

  # 2. The path contains files:
  if ((length(list.files(path)) > 0) != TRUE) {
    stop("but there are no files in it")
  } else{
    cat("and there are files in it\n\n")
  }
}


#' @importFrom RHRV CreateHRVData SetVerbose LoadBeat BuildNIHR FilterNIHR
preparing_analysis <- function(file, rrs, format, easy_options) {
  hrv.data = CreateHRVData()
  hrv.data = SetVerbose(hrv.data, FALSE)

  hrv.data = tryCatch({
    hrv.data = LoadBeat(
      fileType = format,
      HRVData = hrv.data,
      Recordname = file,
      RecordPath = rrs
    )
    hrv.data
  },
  error = function(cond) {
    stop(
      paste(
        "The file \"",
        file,
        "\" could not be loaded. Check if the file is in the correct format; the specified format was \"",
        format,
        "\".",
        sep = ""
      )
    )
  })
  if (easy_options$verbose) {
    message(c("Loading recording ", file))
  }
  hrv.data = BuildNIHR(hrv.data)
  hrv.data = FilterNIHR(hrv.data)
  hrv.data$Beat = hrv.data$Beat[2:nrow(hrv.data$Beat), ]
  hrv.data
}

#Calls an RHRV function with hrv.data after cleaning the parameters
easy_call <- function(hrv.data, mf, ...) {
  args.list = plotrix::clean.args(list(...), mf)
  args.list$HRVData = hrv.data
  do.call(mf, args.list)
}

# Creating time analysis data frames
#' @importFrom RHRV CreateTimeAnalysis
#' @importFrom foreach foreach
#' @importFrom foreach %dopar% %do%
time_analysis <-
  function(format, files, class, rrs2, easy_options, ...) {
    dataFrame = foreach(
      file = files,
      .combine = rbind.data.frame,
      .export = c("preparing_analysis", "easy_call"),
      # .packages = "RHRV",
      .errorhandling = "pass"
    ) %dopar% {
      suppressMessages(library("RHRV", character.only = TRUE))
      hrv.data = preparing_analysis(file = file, rrs = rrs2, format = format,
                                    easy_options = easy_options)
      hrv.data = easy_call(hrv.data, CreateTimeAnalysis, ...)
      results = hrv.data$TimeAnalysis[]
      name_file = list("filename" = file)
      group = list("group" = class)
      # group_name = list("group" = group)
      row_list = c(name_file, results, group)
      as.data.frame(row_list)
    }
    #@todo  ?Remove size column???
    dataFrame
  }

#' @importFrom foreach foreach
#' @importFrom RHRV CreateFreqAnalysis
#' @importFrom RHRV InterpolateNIHR CalculatePSD CalculateEnergyInPSDBands
freq_analysis <-
  function(format, files, class, rrs2, easy_options, ...) {
    dataFrame = foreach(
      file = files,
      .combine = rbind.data.frame,
      .export = c("preparing_analysis", "easy_call"),
      # .packages = "RHRV",
      .errorhandling = "pass"
    ) %dopar% {
      suppressMessages(library("RHRV", character.only = TRUE))
      hrv.data = preparing_analysis(file = file, rrs = rrs2, format = format,
                                    easy_options = easy_options)
      hrv.data = easy_call(hrv.data, InterpolateNIHR, ...)
      zero_indexes = which(hrv.data$HR == 0)
      hr_median = median(hrv.data$HR[-zero_indexes])
      hrv.data$HR[zero_indexes] = hr_median
      hrv.data = easy_call(hrv.data, CreateFreqAnalysis, ...)
      hrv.data = easy_call(hrv.data, CalculatePSD, doPlot = F, ...)
      name_file = list("filename" = file)
      x1 = easy_call(hrv.data, CalculateEnergyInPSDBands, ...)
      names(x1) = c("ULF", "VLF", "LF", "HF")
      group = list("group" = class)
      row_list = c(name_file, x1, group)
      as.data.frame(row_list)
    }
    dataFrame
  }

#' @importFrom foreach foreach
#' @importFrom RHRV CreateFreqAnalysis
#' @importFrom RHRV InterpolateNIHR CalculatePowerBand
wavelet_analysis <-
  function(format, files, class, rrs2, easy_options, ...) {
    dataFrameMWavelet = foreach(
      file = files,
      .combine = rbind.data.frame,
      .export = c("preparing_analysis", "easy_call"),
      # .packages = "RHRV",
      .errorhandling = "pass"
    ) %dopar% {
      suppressMessages(library("RHRV", character.only = TRUE))
      hrv.data = preparing_analysis(file = file, rrs = rrs2, format = format, easy_options = easy_options)
      hrv.data = easy_call(hrv.data, InterpolateNIHR, ...)
      zero_indexes = which(hrv.data$HR == 0)
      hr_median = median(hrv.data$HR[-zero_indexes])
      hrv.data$HR[zero_indexes] = hr_median

      hrv.data = easy_call(hrv.data, CreateFreqAnalysis, ...)
      hrv.data = SetVerbose(hrv.data, easy_options$verbose)
      hrv.data = easy_call(hrv.data, CalculatePowerBand, ...)

      index = length(hrv.data$FreqAnalysis)
      resultsWavelet = hrv.data$FreqAnalysis[[index]]
      resultsWavelet$File = file
      resultsWavelet$HRV = NA
      resultsWavelet$ULF = sum(hrv.data$FreqAnalysis[[index]]$ULF)
      resultsWavelet$VLF = sum(hrv.data$FreqAnalysis[[index]]$VLF)
      resultsWavelet$LF = sum(hrv.data$FreqAnalysis[[index]]$LF)
      resultsWavelet$HF = sum(hrv.data$FreqAnalysis[[index]]$HF)
      resultsWavelet$LFHF = NA
      resultsWavelet$Time = NA
      name_file = list()
      x1 = as.list(resultsWavelet)
      group = list("group" = class)
      row_list = c(name_file, x1, group)
      as.data.frame(row_list)
    }
    dataFrameMWavelet
  }



dunnfreq <- function(dfM, correctionMethod, easy_options = easy_options) {
  dfM$group = factor(dfM$group)
  list(
    ULF = posthoc.kruskal.dunn.test.CheckAllValuesEqual(ULF ~ group, data =
                                                          dfM, p.adjt = correctionMethod, easy_options = easy_options),
    VLF = posthoc.kruskal.dunn.test.CheckAllValuesEqual(VLF ~ group, data =
                                                          dfM, p.adjt = correctionMethod, easy_options = easy_options),
    LF = posthoc.kruskal.dunn.test.CheckAllValuesEqual(LF ~ group, data =
                                                         dfM, p.adjt = correctionMethod, easy_options = easy_options),
    HF = posthoc.kruskal.dunn.test.CheckAllValuesEqual(HF ~ group, data =
                                                         dfM, p.adjt = correctionMethod, easy_options = easy_options)
  )
}

dunntime <- function(dfM, correctionMethod, easy_options = easy_options) {
  dfM$group = factor(dfM$group)
  list(
    SDNN = posthoc.kruskal.dunn.test.CheckAllValuesEqual(SDNN ~ group, data = dfM, p.adjt =
                                                           correctionMethod, easy_options = easy_options),
    SDANN = posthoc.kruskal.dunn.test.CheckAllValuesEqual(SDANN ~ group, data = dfM, p.adjt =
                                                            correctionMethod, easy_options = easy_options),
    SDNNIDX = posthoc.kruskal.dunn.test.CheckAllValuesEqual(SDNNIDX ~ group, data = dfM, p.adjt =
                                                              correctionMethod, easy_options = easy_options),
    pNN50 = posthoc.kruskal.dunn.test.CheckAllValuesEqual(pNN50 ~ group, data = dfM, p.adjt =
                                                            correctionMethod, easy_options = easy_options),
    SDSD = posthoc.kruskal.dunn.test.CheckAllValuesEqual(SDSD ~ group, data = dfM, p.adjt =
                                                           correctionMethod, easy_options = easy_options),
    rMSSD = posthoc.kruskal.dunn.test.CheckAllValuesEqual(rMSSD ~ group, data = dfM, p.adjt =
                                                            correctionMethod, easy_options = easy_options),
    IRRR = posthoc.kruskal.dunn.test.CheckAllValuesEqual(IRRR ~ group, data = dfM, p.adjt =
                                                           correctionMethod, easy_options = easy_options),
    MADRR = posthoc.kruskal.dunn.test.CheckAllValuesEqual(MADRR ~ group, data = dfM, p.adjt =
                                                            correctionMethod, easy_options = easy_options),
    TINN = posthoc.kruskal.dunn.test.CheckAllValuesEqual(TINN ~ group, data = dfM, p.adjt =
                                                           correctionMethod, easy_options = easy_options),
    HRVi = posthoc.kruskal.dunn.test.CheckAllValuesEqual(HRVi ~ group, data = dfM, p.adjt =
                                                           correctionMethod, easy_options = easy_options)
  )
}


#' @importFrom PMCMR posthoc.kruskal.dunn.test
posthoc.kruskal.dunn.test.CheckAllValuesEqual <-
  function(formula, data, p.adjt, easy_options) {
    # If we cannot do the test, I see because all the numerical values are the same,
    # we will return NULL since there are no differences between the populations.
    dunn = tryCatch({
      posthoc.kruskal.dunn.test(formula,
                                data,
                                p.adjust.method =  p.adjt,
                                na.action = na.omit)
    },
    error = function(cond) {
      if (easy_options$verbose) {
        message(c(
          "All values identical in kruskal.dunn.test; pvalue set to 1 for ",
          formula
        ))
      }
      NULL
    })
    dunn
  }

shapiro.test.CheckAllValuesEqual <- function(x, easy_options) {
  # If we cannot do the test, I see because all the numerical values are the same,
  # we will return 1 since there are no differences between the populations.
  pval = tryCatch({
    shapiro.test(x)$p.value
  },
  error = function(cond) {
    message(x)
    if (easy_options$verbose) {
      message(
        "All indice's values identical in shapiro.test, or less than three values are different from NA; non normality assumed and pvalue set to 0"
      )
    }
    0
  })
  pval
}

statistical_analysisFreq <-
  function(dfM,
           numberOfExperimentalGroups,
           correctionMethod,
           easy_options) {
    anova = list(
      ULF = NA,
      VLF = NA,
      LF = NA,
      HF = NA
    )
    kruskal = list(
      ULF = NA,
      VLF = NA,
      LF = NA,
      HF = NA
    )
    dunn = NA
    list = list(anova = anova,
                kruskal = kruskal,
                dunn = dunn)

    listDF = split(dfM, dfM$group)

    dataFramePvalues = data.frame()
    vec = list(
      "group" = NA,
      "p-value ULF" = NA,
      "p-value VLF" = NA,
      "p-value LF" = NA,
      "p-value HF" = NA
    )


    for (objeto in names(listDF)) {
      vec$group = objeto

      for (column in c('ULF', 'VLF', 'LF', 'HF')) {
        destino = paste0('p-value ', column)
        vec[[destino]] = shapiro.test.CheckAllValuesEqual(listDF[[objeto]][[column]], easy_options)
      }

      df = data.frame(vec)

      dataFramePvalues = rbind(dataFramePvalues, df)
    }

    for (column in c('ULF', 'VLF', 'LF', 'HF')) {
      p_values = formula_str = paste0("p.value.", column)
      formula_str = paste0(column, "~ group")
      formula = as.formula(formula_str)

      if (numberOfExperimentalGroups > 2 ||
          all(dataFramePvalues[[p_values]] > easy_options$significance)) {
        if (easy_options$verbose) {
          message(column,
                  " Normal: Anova. P-values = ",
                  dataFramePvalues[[p_values]],
                  "\n")
        }
        list$anova[[column]] = aov(formula, data = dfM, na.action = na.exclude)
      } else {
        if (easy_options$verbose) {
          message(column,
                  " NOT normal: Kruskal. P-values = ",
                  dataFramePvalues[[p_values]],
                  "\n")
        }
        list$kruskal[[column]] = kruskal.test(formula, data = dfM, na.action = na.exclude)
      }

    }

    list$dunn = dunnfreq(dfM, correctionMethod, easy_options = easy_options)
    list

  }

statistical_analysisTime <-
  function(dfM,
           numberOfExperimentalGroups,
           correctionMethod,
           easy_options) {
    anova = list(
      SDNN = NA,
      SDANN = NA,
      SDNNIDX = NA,
      pNN50 = NA,
      SDSD = NA,
      rMSSD = NA,
      IRRR = NA,
      MADRR = NA,
      TINN = NA,
      HRVi = NA
    )
    kruskal = list(
      SDNN = NA,
      SDANN = NA,
      SDNNIDX = NA,
      pNN50 = NA,
      SDSD = NA,
      rMSSD = NA,
      IRRR = NA,
      MADRR = NA,
      TINN = NA,
      HRVi = NA
    )
    dunn = NA
    list = list(anova = anova,
                kruskal = kruskal,
                dunn = dunn)

    listDF = split(dfM, dfM$group)

    dataFramePvalues = data.frame()

    vec = list(
      "group" = NA,
      "p-value SDNN" = NA,
      "p-value SDANN" = NA,
      "p-value SDNNIDX" = NA,
      "p-value pNN50" = NA,
      "p-value SDSD" = NA,
      "p-value rMSSD" = NA,
      "p-value IRRR" = NA,
      "p-value MADRR" = NA,
      "p-value TINN" = NA,
      "p-value HRVi" = NA
    )

    for (objeto in names(listDF)) {
      vec$group = objeto

      for (column in c(
        'SDNN',
        'SDANN',
        'SDNNIDX',
        'pNN50',
        'SDSD',
        'rMSSD',
        'IRRR',
        'MADRR',
        'TINN',
        'HRVi'
      )) {
        destino = paste0('p-value ', column)
        vec[[destino]] =  shapiro.test.CheckAllValuesEqual(listDF[[objeto]][[column]], easy_options)
      }

      df = data.frame(vec)

      dataFramePvalues = rbind(dataFramePvalues, df)
    }

    for (column in c('SDNN',
                     'SDANN',
                     'SDNNIDX',
                     'pNN50',
                     'SDSD',
                     'rMSSD',
                     'IRRR',
                     'MADRR',
                     'TINN',
                     'HRVi')) {
      p_values = formula_str = paste0("p.value.", column)
      formula_str = paste0(column, "~ group")
      formula = as.formula(formula_str)

      if (numberOfExperimentalGroups > 2 ||
          all(dataFramePvalues[[p_values]] > easy_options$significance)) {
        if (easy_options$verbose) {
          cat(column,
              " Normal: Anova. P-values = ",
              dataFramePvalues[[p_values]],
              "\n")
        }
        list$anova[[column]] = aov(formula, data = dfM, na.action = na.exclude)
      } else {
        if (easy_options$verbose) {
          cat(column,
              " NOT normal: Kruskal. P-values = ",
              dataFramePvalues[[p_values]],
              "\n")
        }
        list$kruskal[[column]] = kruskal.test(formula, data = dfM, na.action = na.exclude)
      }
    }

    list$dunn = dunntime(dfM, correctionMethod, easy_options = easy_options)
    list
  }

statistical_analysisNonLinear <-
  function(dfM,
           numberOfExperimentalGroups,
           correctionMethod,
           easy_options) {
    anova = list(
      CorrelationStatistic = NA,
      SampleEntropy = NA,
      MaxLyapunov = NA,
      REC = NA,
      RATIO = NA,
      DET = NA,
      DIV = NA,
      Lmax = NA,
      Lmean = NA,
      LmeanWithoutMain = NA,
      ENTR = NA,
      TREND = NA,
      LAM = NA,
      Vmax = NA,
      Vmean = NA,
      PoincareSD1 = NA,
      PoincareSD2 = NA
    )
    kruskal = list(
      CorrelationStatistic = NA,
      SampleEntropy = NA,
      MaxLyapunov = NA,
      REC = NA,
      RATIO = NA,
      DET = NA,
      DIV = NA,
      Lmax = NA,
      Lmean = NA,
      LmeanWithoutMain = NA,
      ENTR = NA,
      TREND = NA,
      LAM = NA,
      Vmax = NA,
      Vmean = NA,
      PoincareSD1 = NA,
      PoincareSD2 = NA
    )
    dunn = NA
    list = list(anova = anova,
                kruskal = kruskal,
                dunn = dunn)

    listDF = split(dfM, dfM$group)

    dataFramePvalues = data.frame()
    vec = list(
      "group" = NA,
      "p-value CorrelationStatistic" = NA,
      "p-value SampleEntropy" = NA,
      "p-value MaxLyapunov" = NA,
      "p-value REC" = NA,
      "p-value RATIO" = NA,
      "p-value DET" = NA,
      "p-value DIV" = NA,
      "p-value Lmax" = NA,
      "p-value Lmean" = NA,
      "p-value LmeanWithoutMain" = NA,
      "p-value ENTR" = NA,
      "p-value TREND" = NA,
      "p-value LAM" = NA,
      "p-value Vmax" = NA,
      "p-value Vmean" = NA,
      "p-value PoincareSD1" = NA,
      "p-value PoincareSD2" = NA
    )

    for (objeto in names(listDF)) {
      vec$group = objeto

      for (column in c(
        'CorrelationStatistic',
        'SampleEntropy',
        'MaxLyapunov',
        'REC',
        'RATIO',
        'DET',
        'DIV',
        'Lmax',
        'Lmean',
        'LmeanWithoutMain',
        'ENTR',
        'TREND',
        'LAM',
        'Vmax',
        'Vmean',
        'PoincareSD1',
        'PoincareSD2'
      )) {
        destino = paste0('p-value ', column)
        vec[[destino]] =  shapiro.test.CheckAllValuesEqual(listDF[[objeto]][[column]], easy_options)
      }

      df = data.frame(vec)

      dataFramePvalues = rbind(dataFramePvalues, df)
    }

    for (column in c(
      'CorrelationStatistic',
      'SampleEntropy',
      'MaxLyapunov',
      'REC',
      'RATIO',
      'DET',
      'DIV',
      'Lmax',
      'Lmean',
      'LmeanWithoutMain',
      'ENTR',
      'TREND',
      'LAM',
      'Vmax',
      'Vmean',
      'PoincareSD1',
      'PoincareSD2'
    )) {
      p_values = formula_str = paste0("p.value.", column)
      formula_str = paste0(column, "~ group")
      formula = as.formula(formula_str)
      if (numberOfExperimentalGroups > 2 ||
          all(dataFramePvalues[[p_values]] > easy_options$significance)) {
        if (easy_options$verbose) {
          cat(column,
              " Normal: Anova. P-values = ",
              dataFramePvalues[[p_values]],
              "\n")
        }

        #ANOVA will fail if the hrv statistic could not be calculated for all recordings in a group
        list$anova[[column]] = tryCatch({
          aov(formula, data = dfM, na.action = na.exclude)
        },
        error = function(cond) {
          NA
        })


      } else {
        if (easy_options$verbose) {
          cat(column,
              " NOT normal: Kruskal. P-values = ",
              dataFramePvalues[[p_values]],
              "\n")
        }
        #Krustal will fail if the statistic could not be calculated for all recordings in a group
        list$kruskal[[column]] = tryCatch({
          kruskal.test(formula, data = dfM, na.action = na.exclude)
        },
        error = function(cond) {
          NA
        })
      }
    }
    list$dunn = dunnNonLinear(dfM, correctionMethod, easy_options = easy_options)
    list

  }


colectpValues <-
  function(listTime, listFreq, listNonLinear) {
    listpValues = list(
      ULF = NA,
      VLF = NA,
      LF = NA,
      HF = NA,
      SDNN = NA,
      SDANN = NA,
      SDNNIDX = NA,
      pNN50 = NA,
      SDSD = NA,
      rMSSD = NA,
      IRRR = NA,
      MADRR = NA,
      TINN = NA,
      HRVi = NA,
      CorrelationStatistic = NA,
      SampleEntropy = NA,
      MaxLyapunov = NA,
      REC = NA,
      RATIO = NA,
      DET = NA,
      DIV = NA,
      Lmax = NA,
      Lmean = NA,
      LmeanWithoutMain = NA,
      ENTR = NA,
      TREND = NA,
      LAM = NA,
      Vmax = NA,
      Vmean = NA,
      PoincareSD1 = NA,
      PoincareSD2 = NA
    )

    for (column in c('ULF', 'VLF', 'LF', 'HF')) {
      if (!inherits(listFreq$anova[[column]], "aov")) {
        listpValues[[column]] = listFreq$kruskal[[column]]$p.value
      } else{
        listpValues[[column]] = extract_ANOVA_pvalue(listFreq$anova[[column]])
      }
    }

    for (column in c('SDNN',
                     'SDANN',
                     'SDNNIDX',
                     'pNN50',
                     'SDSD',
                     'rMSSD',
                     'IRRR',
                     'MADRR',
                     'TINN',
                     'HRVi')) {
      if (!inherits(listTime$anova[[column]], "aov")) {
        listpValues[[column]] = listTime$kruskal[[column]]$p.value
      } else{
        listpValues[[column]] = extract_ANOVA_pvalue(listTime$anova[[column]])
      }
    }

    # In order for it to only be performed when there is non linear results:
    if (!all(is.na(listNonLinear))) {
      for (column in c(
        'CorrelationStatistic',
        'SampleEntropy',
        'MaxLyapunov',
        'REC',
        'RATIO',
        'DET',
        'DIV',
        'Lmax',
        'Lmean',
        'LmeanWithoutMain',
        'ENTR',
        'TREND',
        'LAM',
        'Vmax',
        'Vmean',
        'PoincareSD1',
        'PoincareSD2'
      )) {
        if (is.na(listNonLinear[["anova"]][[column]])) {
          if (inherits(listNonLinear[["kruskal"]][[column]], "htest")) {
            listpValues[[column]] = listNonLinear[["kruskal"]][[column]]$p.value
          }
          else{
            #if we have not been able to calculate the statistic,
            #we cannot affirm that there are differences in the statistic
            listpValues[[column]] = 1
          }
        } else{
          p.val.tmp = extract_ANOVA_pvalue(listNonLinear[["anova"]][[column]])
          if (is.na(p.val.tmp)) {
            #if we have not been able to calculate the statistic,
            #we cannot affirm that there are differences in the statistic
            listpValues[[column]] = 1
          }
          else{
            listpValues[[column]] = extract_ANOVA_pvalue(listNonLinear[["anova"]][[column]])
          }
        }
      }
    }
    listpValues
  }

correctpValues <-
  function(listpValues,
           correction,
           correctionMethod) {
    listpValuesCorrected = list(
      ULF = NA,
      VLF = NA,
      LF = NA,
      HF = NA,
      SDNN = NA,
      SDANN = NA,
      SDNNIDX = NA,
      pNN50 = NA,
      SDSD = NA,
      rMSSD = NA,
      IRRR = NA,
      MADRR = NA,
      TINN = NA,
      HRVi = NA,
      CorrelationStatistic = NA,
      SampleEntropy = NA,
      MaxLyapunov = NA,
      REC = NA,
      RATIO = NA,
      DET = NA,
      DIV = NA,
      Lmax = NA,
      Lmean = NA,
      LmeanWithoutMain = NA,
      ENTR = NA,
      TREND = NA,
      LAM = NA,
      Vmax = NA,
      Vmean = NA,
      PoincareSD1 = NA,
      PoincareSD2 = NA
    )

    if (correction == TRUE) {
      listpValuesCorrected = p.adjust(listpValues, correctionMethod)
      listpValuesCorrected <- as.list(listpValuesCorrected)

    } else{
      listpValuesCorrected = listpValues
    }
    listpValuesCorrected
  }

split_path <- function(path) {
  if (dirname(path) %in% c(".", path))
    return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

extract_ANOVA_pvalue <- function(anovaObject) {
  pvalue = summary(anovaObject)[[1]][1, 5]
  pvalue
}

#' @export
print.RHRVEasyResult <- function(results) {
  easy_options <- attr(results, "easy_options")
  listDF = split(results$TimeAnalysis, results$TimeAnalysis$group)

  differencesFound = FALSE

  cat(
    "\n\nResult of the analysis of the variability of the heart rate of the group",
    levels(results$TimeAnalysis$group)[1],
    "versus the group",
    levels(results$TimeAnalysis$group)[2],
    ":\n\n"
  )

  for (column in c('SDNN',
                   'SDANN',
                   'SDNNIDX',
                   'pNN50',
                   'SDSD',
                   'rMSSD',
                   'IRRR',
                   'MADRR',
                   'TINN',
                   'HRVi')) {
    if (all(is.na(results$StatysticalAnalysisTime$anova[[column]]))) {
      #report kruskal
      if (!is.na(results$pValues[[column]]) &&
          results$pValues[[column]] < easy_options$significance) {
        #error pvalue 1
        differencesFound = TRUE
        cat(
          "\nThere is a statistically significant difference in",
          column,
          "; pvalue: ",
          results$pValues[[column]],
          "\n"
        )

        for (i in 1:length(listDF)) {
          group = names(listDF)[i]
          cat(
            column,
            " for the group",
            levels(results$TimeAnalysis$group)[i],
            "is",
            mean(listDF[[group]][[column]], na.rm = TRUE),
            "+-",
            sd(listDF[[group]][[column]], na.rm = TRUE),
            "\n"
          )
        }
      }
    }
    #report anova
    else{
      if (!is.na(results$pValues[[column]]) &&
          results$pValues[[column]] < easy_options$significance) {
        differencesFound = TRUE
        cat(
          "\nThere is a statistically significant difference in",
          column,
          "; pvalue: ",
          results$pValues[[column]],
          "\n"
        )

        for (i in 1:length(listDF)) {
          group = names(listDF)[i]
          cat(
            column,
            " for the group ",
            levels(results$TimeAnalysis$group)[i],
            "is",
            mean(listDF[[group]][[column]], na.rm = TRUE),
            "+-",
            sd(listDF[[group]][[column]], na.rm = TRUE),
            "\n"
          )
        }

      }
    }


    # Two Conditions to report Dunn:
    # 1. We have more than 2 groups, we check that by looking at the length of listDF
    if (length(listDF) > 2) {
      # 2. ANOVA Test is significative. We check that by comparing it to the signif_level
      var = which(results$StatysticalAnalysisTime$dunn[[column]][["p.value"]] <
                    easy_options$significance,
                  arr.ind = TRUE)

      if (length(var) > 0) {
        cat(
          "\nGroups with stastically significant differences in ",
          column,
          " according to the Dunn test :\n"
        )
        print(results$StatysticalAnalysisTime$dunn[[column]])
      }
    }

  }

  listDF = split(results$FrequencyAnalysis,
                 results$FrequencyAnalysis$group)

  for (column in c('ULF', 'VLF', 'LF', 'HF')) {
    if (all(is.na(results$StatysticalAnalysisFrequency$anova[[column]]))) {
      #report kruskal
      if (results$pValues[[column]] < easy_options$significance) {
        #error pvalue 1
        differencesFound = TRUE
        cat(
          "\nThere is a statistically significant difference in",
          column,
          "; pvalue: ",
          results$pValues[[column]],
          "\n"
        )

        for (i in 1:length(listDF)) {
          group = names(listDF)[i]
          cat(
            column,
            " for the group",
            levels(results$TimeAnalysis$group)[i],
            "is",
            mean(listDF[[group]][[column]], na.rm = TRUE),
            "+-",
            sd(listDF[[group]][[column]], na.rm = TRUE),
            "\n"
          )
        }
      }
    }
    #report anova
    else{
      if (results$pValues[[column]] < easy_options$significance) {
        differencesFound = TRUE
        cat(
          "\nThere is a statistically significant difference in",
          column,
          "; pvalue: ",
          results$pValues[[column]],
          "\n"
        )

        for (i in 1:length(listDF)) {
          group = names(listDF)[i]
          cat(
            column,
            " for the group ",
            levels(results$TimeAnalysis$group)[i],
            "is",
            mean(listDF[[group]][[column]], na.rm = TRUE),
            "+-",
            sd(listDF[[group]][[column]], na.rm = TRUE),
            "\n"
          )
        }

      }
    }

    # Two Conditions to report Dunn:
    # 1. We have more than 2 groups, we check that by looking at the length of listDF

    if(length(listDF)>2){

      #@TODO creo que deber?amos cambiar la condici?n para reportar Dunn por
      #  if(results$pValues[[column]]<signif_level). Si no en alguna ocasi?n (LF en la siguiente prueba)
      #reporta Dunn sin haber reportado ANOVA
      #a2b=RHRVEasy(folders =c("C:\\rrs\\RHRVEasy\\rrs\\normal",
      #                        "C:\\rrs\\RHRVEasy\\rrs\\chf",
      #                        "C:\\rrs\\RHRVEasy\\rrs\\normal_half",
      #                        "C:\\rrs\\RHRVEasy\\rrs\\chf_half"), significance_level = 0.05)

      #@TODO Lo mismo en otros sitios
      # 2. ANOVA Test is significative. We check that by comparing it to the signif_level

      variable = which(results$StatysticalAnalysisFrequency$dunn[[column]][["p.value"]] <
                         easy_options$significance,
                       arr.ind = TRUE)
      if (length(variable) > 0) {
        cat(
          "\nGroups with stastically significant differences in ",
          column,
          " according to the Dunn test :\n"
        )
        print(results$StatysticalAnalysisFrequency$dunn[[column]])
      }

    }

  }

  if (!all(is.na(results$NonLinearAnalysis))) {
    listDF = split(results$NonLinearAnalysis,
                   results$NonLinearAnalysis$group)

    for (column in c(
      'CorrelationStatistic',
      'SampleEntropy',
      'MaxLyapunov',
      'REC',
      'RATIO',
      'DET',
      'DIV',
      'Lmax',
      'Lmean',
      'LmeanWithoutMain',
      'ENTR',
      'TREND',
      'LAM',
      'Vmax',
      'Vmean',
      'PoincareSD1',
      'PoincareSD2'
    )) {
      if (all(is.na(results$StatysticalAnalysisNonLinear$anova[[column]]))) {
        #report kruskal
        if (results$pValues[[column]] < easy_options$significance) {
          #error pvalue 1
          differencesFound = TRUE
          cat(
            "\nThere is a statistically significant difference in",
            column,
            "; pvalue: ",
            results$pValues[[column]],
            "\n"
          )

          for (i in 1:length(listDF)) {
            group = names(listDF)[i]
            cat(
              column,
              " for the group",
              levels(results$TimeAnalysis$group)[i],
              "is",
              mean(listDF[[group]][[column]], na.rm = TRUE),
              "+-",
              sd(listDF[[group]][[column]], na.rm = TRUE),
              "\n"
            )
          }
        }
      }
      #report anova
      else{
        if (results$pValues[[column]] < easy_options$significance) {
          differencesFound = TRUE
          cat(
            "\nThere is a statistically significant difference in",
            column,
            "; pvalue: ",
            results$pValues[[column]],
            "\n"
          )

          for (i in 1:length(listDF)) {
            group = names(listDF)[i]
            cat(
              column,
              " for the group ",
              levels(results$TimeAnalysis$group)[i],
              "is",
              mean(listDF[[group]][[column]], na.rm = TRUE),
              "+-",
              sd(listDF[[group]][[column]], na.rm = TRUE),
              "\n"
            )
          }

        }
      }

      # Two Conditions to report Dunn:
      # 1. We have more than 2 groups, we check that by looking at the length of listDF

      if (length(listDF) > 2) {
        # 2. ANOVA Test is significative. We check that by comparing it to the signif_level

        variable = which(results$StatysticalAnalysisNonLinear$dunn[[column]][["p.value"]] <
                           easy_options$significance,
                         arr.ind = TRUE)
        if (length(variable) > 0) {
          cat(
            "\nGroups with stastically significant differences in ",
            column,
            " according to the Dunn test :\n"
          )
          print(results$StatysticalAnalysisNonLinear$dunn[[column]])
        }

      }

    }
  }

  if (!differencesFound) {
    cat("No statistically significant difference were found\n")
  }
}

#' @importFrom writexl write_xlsx
#' @export
saveHRVindexes <- function(results, saveHRVindexesInPath = ".") {
  #if called directly by RHRVEasy witout an explicit value for saveHRVindexesInPath
  #then saveHRVindexesInPath is null ans nothing hapens
  if (!is.null(saveHRVindexesInPath)) {
    tryCatch({
      me = merge(results$TimeAnalysis, results$FrequencyAnalysis)
      if (ncol(results$NonLinearAnalysis) == 0) {
        #no linear analysis
        frameTosave = me
      } else{
        frameTosave = merge(me, results$NonLinearAnalysis)
      }
      fileName = ""
      for (lev in levels(as.factor(results$TimeAnalysis$group))) {
        fileName = paste(fileName, lev, " vs ", sep = "")
      }
      fileName = substr(fileName, 1, nchar(fileName) - 4)

      fileName = paste(saveHRVindexesInPath, "/", fileName, ".xlsx", sep =
                         "")
      write_xlsx(frameTosave, fileName)
    },
    error = function(cond) {
      message("There was an error when trying to save the results to the Ecel file")
      message(cond)
    })
  }
}


#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
.prepare_cluster <- function(n_jobs, verbose) {
  n_cores <- parallel::detectCores(logical = FALSE)
  if (n_jobs <= 0) {
    n_jobs <- n_cores
  } else if (n_jobs > n_cores) {
    n_jobs <- n_cores
  }
  if (n_jobs > 1) {
    cl <- parallel::makeCluster(n_jobs, outfile="") # using outfile = "" may be useful for debugging
    doParallel::registerDoParallel(cl)
    if (verbose) {
      message(paste("Registering cluster with", n_jobs, "nodes"))
    }
  } else {
    cl <- NULL
    foreach::registerDoSEQ()
  }
  cl
}


.build_easy_options <- function(verbose, significance) {
  list(
    "verbose" = verbose,
    "significance" = significance
  )
}


#' @export
RHRVEasy <-
  function(folders,
           correction = FALSE,
           correctionMethod = "bonferroni",
           verbose = FALSE,
           format = "RR",
           typeAnalysis = 'fourier',
           significance_level = 0.25,
           nonLinear = FALSE,
           doRQA = FALSE,   # ignored if nonLinear = FALSE
           saveHRVindexesInPath = NULL,
           n_jobs = 1,
           ...) {
    dataFrameMWavelet = data.frame()
    dataFrameMTime = data.frame()
    dataFrameMFreq = data.frame()
    dataFrameMNonLinear = data.frame()
    listNonLinearStatisticalAnalysis = list()
    listTimeStatysticalAnalysis = list()
    listFreqStatysticalAnalysis = list()

    cl <- .prepare_cluster(n_jobs, verbose)

    easy_options <- .build_easy_options(
      verbose = verbose,
      significance = significance_level
    )

    for (folder in folders) {
      file_validation(folder)
      dataFrameMTime = rbind(
        dataFrameMTime,
        dataFrameMTime = time_analysis(
          format,
          list.files(folder),
          split_path(folder)[1],
          folder,
          easy_options,
          ...
        )
      )
      if (nonLinear) {
        if (easy_options$verbose) {
          message("Performing non linear analysis...")
        }
        dataFrameMNonLinear = rbind(
          dataFrameMNonLinear,
          non_linear_analysis(
            format,
            list.files(folder),
            split_path(folder)[1],
            folder,
            easy_options,
            doRQA,
            ...
          )
        )
      }
    }

    numberOfExperimentalGroups = length(folders)
    # Statistical analysis of both

    listTimeStatysticalAnalysis = statistical_analysisTime(dataFrameMTime,
                                                           numberOfExperimentalGroups,
                                                           correctionMethod,
                                                           easy_options)

    # FREQUENCY:
    if (typeAnalysis == "fourier") {
      for (folder in folders) {
        dataFrameMFreq = rbind(
          dataFrameMFreq,
          dataFrameMFreq = freq_analysis(
            format,
            list.files(folder),
            split_path(folder)[1],
            folder,
            easy_options,
            ...
          )
        )
      }

      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMFreq,
                                                             numberOfExperimentalGroups,
                                                             correctionMethod,
                                                             easy_options)
    }

    # WAVELET
    if (typeAnalysis == "wavelet") {
      for (folder in folders) {
        dataFrameMWavelet = rbind(
          dataFrameMWavelet,
          dataFrameMWavelet =
            wavelet_analysis(
              format,
              list.files(folder),
              split_path(folder)[1],
              folder,
              type = typeAnalysis,
              easy_options,
              ...
            )
        )
      }
      listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMWavelet,
                                                             numberOfExperimentalGroups,
                                                             correctionMethod,
                                                             easy_options)

      dataFrameMFreq = dataFrameMWavelet
    }

    if (!is.null(cl)) {
      parallel::stopCluster(cl)
    }

    if (!all(is.na(dataFrameMNonLinear))) {
      listNonLinearStatisticalAnalysis = statistical_analysisNonLinear(
        dataFrameMNonLinear,
        numberOfExperimentalGroups,
        correctionMethod,
        easy_options
      )
    } else{
      listNonLinearStatisticalAnalysis = NA
    }

    # FIXME
    uncorrectedPvalues = colectpValues(
      listTimeStatysticalAnalysis,
      listFreqStatysticalAnalysis,
      listNonLinearStatisticalAnalysis
    )
    listpValues = correctpValues(uncorrectedPvalues, correction, correctionMethod)

    results = list(
      "TimeAnalysis" = dataFrameMTime,
      "StatysticalAnalysisTime" = listTimeStatysticalAnalysis,
      "FrequencyAnalysis" = dataFrameMFreq,
      "StatysticalAnalysisFrequency" = listFreqStatysticalAnalysis,
      "NonLinearAnalysis" = dataFrameMNonLinear,
      "StatysticalAnalysisNonLinear" = listNonLinearStatisticalAnalysis,
      "pValues" = listpValues,
      "uncorrectedPvalues" = uncorrectedPvalues
    )
    class(results) = "RHRVEasyResult"
    attr(results, "easy_options") <- easy_options
    saveHRVindexes(results, saveHRVindexesInPath)
    results

  }