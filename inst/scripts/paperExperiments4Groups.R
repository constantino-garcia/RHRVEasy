devtools::load_all(".")

# easyAnalysis <- RHRVEasy(
#  folders = c("rrs/chf/", "rrs/normal/", "rrs/chf_half/", "rrs/normal_half/"),
#  nonLinear = TRUE,
#  doRQA = TRUE,
#  verbose = TRUE,
#  nJobs = 25,
#  typeAnalysis = "fourier"
# )
#saveRDS(easyAnalysis, "paperExperiments.RDS")

easyAnalysis <- readRDS("paperExperiments.RDS")

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

saveRDS(redoStats(easyAnalysis, "fdr"), "paperExperimentsFDR.RDS")
saveRDS(redoStats(easyAnalysis, "none"), "paperExperimentsNONE.RDS")


library("tidyverse")
 easyAnalysis <- readRDS("paperExperiments.RDS")
 bonferroni <- readRDS("paperExperiments.RDS")
 fdr <- readRDS("paperExperimentsFDR.RDS")
 none <- readRDS("paperExperimentsNONE.RDS")

 stopifnot(all(bonferroni$HRVIndices == fdr$HRVIndices, na.rm = TRUE))
 stopifnot(all(none$HRVIndices == fdr$HRVIndices, na.rm = TRUE))

 bonferroni <- dplyr::rename(bonferroni$stats, bonferroni = adj.p.value)
 fdr <- dplyr::rename(fdr$stats, fdr = adj.p.value)
 none <- dplyr::rename(none$stats, none = adj.p.value)

 pvalues <- merge(dplyr::select(bonferroni, -pairwise), dplyr::select(fdr, -pairwise))
 stopifnot(all(none$p.value == none$none))
 none <- dplyr::arrange(none, HRVIndex)
 pvalues <- dplyr::arrange(pvalues, HRVIndex)
 stopifnot(all(none$p.value == pvalues$p.value))

 significance <- 0.05
 pvalues %>% filter(p.value < significance, bonferroni < significance, fdr < significance) %>% pull("HRVIndex")
