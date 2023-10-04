#debugSource("R/nonlinearAnalysis.R")
#source("R/ScalingRegionEstimation.R")
#source("R/RHRVEasy.R")
library("RHRVEasy")

# TODO devtools::install_github("constantino-garcia/nonlinearTseries", ref = "rqa")
#library("RHRV")
#library("foreach")



results <- RHRVEasy(
  folders = c("/data/research/rhrveasy/current_chf", "/data/research/rhrveasy/current_normal/"),
  nonLinear = TRUE,
  doRQA = FALSE,
  verbose = TRUE,
  n_jobs = 2
)

print(results)
