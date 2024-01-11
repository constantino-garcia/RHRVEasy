devtools::load_all(".")
#library("RHRVEasy")
rm(anal)
anal <- RHRVEasy(
  folders = c("/data/easy/current_chf2/", "/data/easy/current_normal2/"),
  nonLinear = TRUE,
  doRQA = TRUE,
  verbose = TRUE,
  nJobs = 2,
  typeAnalysis = "wavelet"
)

print(anal$HRVIndices)

