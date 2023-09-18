devtools::load_all(".")
library("RHRV")


results <- RHRVEasy(
  folders = c("/data/easy/chf_micro2", "/data/easy/normal_micro/"),
  nonLinear = TRUE,
  doRQA = TRUE,
  verbose = TRUE,
  #n_jobs = 1
)

