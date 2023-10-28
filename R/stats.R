tryShapiroTest <- function(x) {
  tryCatch({
    shapiro.test(x)$p.value
  }, error = function(e) {
    # shapiro may fail if we have few data or all indices are equal.
    # Let's assume normality returning a pval of 1
    1
   })
}

#' @importFrom broom tidy
columnStats <- function(HRVIndices, column, easy_options) {
  shapiroPVals = sapply(
    split(results[[column]], results$group),
    tryShapiroTest
  )

  analysisFormula = formula(paste(column, "~ group"))
  if (any(shapiroPVals < easy_options$significance)) {
    # Non-normal data
    statsResults = kruskal.test(analysisFormula, data = HRVIndices)
    statsResults = tidy(statsResults)
  } else {
    # Normal data
    statsResults = aov(analysisFormula, data = HRVIndices)
    statsResults = tidy(statsResults)
    statsResults = statsResults[statsResults$term == "group", ]
    statsResults$method = "ANOVA"
  }
  statsResults = statsResults[, c("p.value", "method")]
  statsResults$HRVIndex = column
  statsResults
}

statsAnalysis <- function(HRVIndices, easy_options) {
  indicesNames = setdiff(colnames(HRVIndices), c("file", "group"))
  pVals = foreach(
    column = indicesNames,
    .combine = rbind
  ) %do% {
    columnStats(HRVIndices, column, easy_options)
  }
  if (!is.null(easy_options$method)) {
    pVals$adj.p.val = p.adjust(pVals$p.value, method = easy_options$method)
  } else {
    pVals$adj.p.val = pVals$p.value
  }
  numberGroups = length(unique(HRVIndices$group))
  # Do post-hoc if numberGroups > 2
  if (numberGroups > 2) {
  } else {
  }
}
