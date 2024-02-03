devtools::load_all(".")


# without nonlinear -------------------------------------------------------

eato <- RHRVEasy(
        folders = c("RRData/chf/", "RRData/normal/"),
        verbose = FALSE,
        nJobs = 25,
        typeAnalysis = "fourier"
)

eatoFdr <- RHRVEasyStats(eato, "fdr")
eatoFdr


mergePVals <- function(bonf, fdr) {
  merge(
    bonf$stats,
    fdr$stats,
    by = setdiff(colnames(bonf$stats), "adj.p.value"),
    suffixes = c(".bonf", ".fdr")
  )
}
print(
  mergePVals(eato, eatoFdr) %>%
    filter(adj.p.value.fdr < 0.05) %>%
    mutate(bonfS = (adj.p.value.bonf < 0.05))
)

# easyAnalysis <- RHRVEasy(
#         folders = c("rrs/chf/", "rrs/normal/"),
#         nonLinear = TRUE,
#         doRQA = TRUE,
#         verbose = TRUE,
#         nJobs = 25,
#         typeAnalysis = "fourier"
# )
# saveRDS(easyAnalysis, "2-paperExperiments.RDS")
#

easyAnalysis <- readRDS("tmp/2-paperExperiments.RDS")

RHRVEasyStats(easyAnalysis, "fdr")
