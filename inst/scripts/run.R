devtools::load_all(".")
# library("RHRVEasy")

# debugonce(nlaSingleFile)
debugonce(RHRVEasy)
anal <- RHRVEasy(
  folders = c("/data/easy/current_chf2/", "/data/easy/current_normal2/"),
  nonLinear = TRUE,
  doRQA = FALSE,
  verbose = TRUE,
  nJobs = 1
)

anal
gg
