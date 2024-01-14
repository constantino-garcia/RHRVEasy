devtools::load_all(".")

easyAnalysis <- RHRVEasy(
        folders = c("rrs/chf/", "rrs/normal/"),
        nonLinear = TRUE,
        doRQA = TRUE,
        verbose = TRUE,
        nJobs = 25,
        typeAnalysis = "fourier"
)
saveRDS(easyAnalysis, "2-paperExperiments.RDS")

easyAnalysis <- readRDS("2-paperExperiments.RDS")

# TODO: move to RHRVEasy
redoStats <- function(easyAnalysis, method) {
        HRVIndices <- easyAnalysis$HRVIndices
        easyOptions <- attr(easyAnalysis, "easyOptions")
        easyOptions$method <- method
        pVals <- statsAnalysis(HRVIndices, easyOptions)
        results <- list("HRVIndices" = HRVIndices, "stats" = pVals)
        class(results) <- "RHRVEasyResult"
        attr(results, "easyOptions") <- easyOptions
        results
}

saveRDS(redoStats(easyAnalysis, "fdr"), "2-paperExperimentsFDR.RDS")
saveRDS(redoStats(easyAnalysis, "none"), "2-paperExperimentsNONE.RDS")
