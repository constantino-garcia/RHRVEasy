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

cd = hrv.data$NonLinearAnalysis[[1]]$correlation
plot(cd$computations)


cd_min = 0.0303
debugonce(RQA)
hrv.data = RQA(
      hrv.data,
      indexNonLinearAnalysis = 1,
      embeddingDim = 4,
      timeLag = 2,
      radius = 10,
      doPlot = FALSE
)



