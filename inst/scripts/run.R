# devtools::load_all(".")
library("RHRVEasy")
rm(anal)
anal <- RHRVEasy(
  folders = c("/data/easy/chf_mini/", "/data/easy/normal_mini//"),
  correctionMethod = "fdr",
  nonLinear = FALSE,
  doRQA = TRUE,
  verbose = TRUE,
  nJobs = 2,
  typeAnalysis = "wavelet",
  saveHRVIndicesInPath = NULL
)

print(anal)

print(anal$HRVIndices)
anal$stats$adj.p.value <- anal$stats$adj.p.value/20

