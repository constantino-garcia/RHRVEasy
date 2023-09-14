devtools::load_all(".")
library("RHRV")

# results <- RHRVEasy(
#   folders = c("/data/easy/chf_micro", "/data/easy/normal_mini/"),
#   nonLinear = TRUE,
#   doRQA = TRUE,
#   verbose = TRUE,
#   #n_jobs = 1
# )

hrv.data = CreateHRVData()
hrv.data = LoadBeatRR(hrv.data, "/data/easy/chf_micro/chf205_rr_secs.txt")
hrv.data = BuildNIHR(hrv.data)
hrv.data = CreateNonLinearAnalysis(hrv.data)
hrv.data = FilterNIHR(hrv.data)
PlotNIHR(hrv.data)

kTimeLag = 2
kTheilerWindow = 100

if (FALSE) {
  hrv.data = CalculateCorrDim(
    hrv.data,
    indexNonLinearAnalysis = 1,
    minEmbeddingDim = 2,
    maxEmbeddingDim = 5,
    timeLag = 2,
    minRadius = 10,
    maxRadius = 50,
    pointsRadius = 20,
    theilerWindow = 10,
    corrOrder = 2,
    doPlot = FALSE
  )
  saveRDS(hrv.data, "hrvWithCD.RDS")
} else {
  hrv.data = readRDS("hrvWithCD.RDS")
}

cd = hrv.data$NonLinearAnalysis[[1]]$correlation$computations

rqaEmbedding = 5
refCorrelation = 1e-3

correlations = cd$corr.matrix[cd$embedding.dims == rqaEmbedding,]
interpolatedFun = approxfun(cd$radius, correlations - refCorrelation)
smallCorrRadius = tryCatch(
  uniroot(
    interpolatedFun, interval = range(cd$radius)
  ),
  error = function(e) e
)
if (inherits(smallCorrRadius, "error")) {
  # warning(paste("Could not find a small radius for neighbor search in file", file))
  smallCorrRadius = min(cd$radius)
  refCorrelation = correlations[cd$radius == smallCorrRadius]
} else {
  smallCorrRadius = smallCorrRadius$root
}
print(paste(" small radius:", smallCorrRadius))
nTakens = length(hrv.data$Beat$RR) - (rqaEmbedding - 1) * kTimeLag
# 4B per entry
estimatedRQASize = (refCorrelation * nTakens * nTakens) * 4

print(estimatedRQASize * 1e-9)
# TODO:Test on windows, mac, etc
availableRam = as.numeric(benchmarkme::get_ram())
if (estimatedRQASize > 0.5 * availableRam) {
  print("Nope")
}

# cd_min = 0.0303
# debugonce(RQA)
# hrv.data = RQA(
#       hrv.data,
#       indexNonLinearAnalysis = 1,
#       embeddingDim = 4,
#       timeLag = 2,
#       radius = 10,
#       doPlot = FALSE
# )
#
#
#
