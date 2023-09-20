devtools::load_all(".")
library("RHRV")


results <- RHRVEasy(
  folders = c("/data/easy/chf", "/data/easy/normal"),
  nonLinear = TRUE,
  doRQA = TRUE,
  verbose = TRUE,
  #n_jobs = 1
)

beepr::beep(8)
