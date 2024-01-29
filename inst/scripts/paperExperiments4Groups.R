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
 HRVIndices <- bonferroni$HRVIndices
 bonferroni <- bonferroni$stats
 fdr <- fdr$stats
 none <- none$stats

 stopifnot(all(none$p.value == none$adj.p.value))
 stopifnot(all(none$p.value == bonferroni$p.value))
 stopifnot(all(none$p.value == fdr$p.value))

getSignificantPosthocs <- function(df, significance = 0.05) {
       posthoc <- do.call(rbind, df[df$adj.p.value < significance, ]$pairwise)
       posthoc <- posthoc[posthoc$adj.p.value < significance, ]
       posthoc$adj.method <- deparse(substitute(df))
       posthoc
}

filterByComparison <- function(posthoc, g1, g2) {
        posthoc %>% filter(((group1 == g1) & (group2 == g2)) | ((group1 == g2) & (group2 == g1)))
}

sb <- getSignificantPosthocs(bonferroni)
sf <- getSignificantPosthocs(fdr)
sn <- getSignificantPosthocs(none)


fdr_sign <- sf %>% filterByComparison("normal", "chf") %>% pull("HRVIndex")
none_sign <- sn %>% filterByComparison("normal", "chf") %>% pull("HRVIndex")
always <- intersect(
        intersect(
                sb %>% filterByComparison("normal", "chf") %>% pull("HRVIndex"),
                fdr_sign
        ),
        none_sign
)

diffIndex <- setdiff(fdr_sign, always)
diffIndex2 <- setdiff(none_sign, always)

sf %>% dplyr::filter(HRVIndex %in% diffIndex) %>%
        filter(group1 == "normal", group2 == "chf")
bonferroni %>% filter(HRVIndex %in% diffIndex)
filter(
        (bonferroni %>% filter(HRVIndex == "SDNNIDX") %>% pull("pairwise"))[[1]],
        group1 == "normal", group2 == "chf"
)

sn %>% dplyr::filter(HRVIndex %in% diffIndex2) %>%
        filter(group1 == "normal", group2 == "chf")


# -------------------------------------------------------------------------
getSignificantGroups <- function(significantDataframe) {
        significantDataframe %>%
                mutate(comparison = paste(sep = " Vs ", group1, group2)) %>%
                group_by(comparison, adj.method) %>%
                summarise(numberSignificantIndices = n())
}
significantGroups <- getSignificantGroups(rbind(sn, sf, sb)) %>%
        pivot_wider(names_from = adj.method, values_from = numberSignificantIndices)
significantGroups <-
        significantGroups %>%
        relocate(none, .before = bonferroni) %>%
        relocate(fdr, .before = bonferroni)
print("Number of significant indices")
print(significantGroups)
knitr::kable(significantGroups, "latex")

significantGroups %>% rowwise() %>% mutate(m = median(fdr, bonferroni))
significantGroups %>% ungroup() %>%  summarise(across(where(is.numeric), mean))

significant_indices <- function(df) {
        d <- df %>%
                mutate(comparison = paste(sep = " Vs ", group1, group2)) %>%
                group_by(comparison, adj.method) %>%
                select(comparison, HRVIndex) %>% nest(data=HRVIndex)
        referenceIndices <- filter(d, comparison == "normal Vs chf") %>% pull("data")
        referenceIndices <- referenceIndices[[1]]$HRVIndex
        for (i in 1:4) {
                        if(length(setdiff(d$data[[i]]$HRVIndex, referenceIndices)) != 0) {
                               stop(d$comparison[[i]], "\n")
                        }
                }
        }
significant_indices(sn)
significant_indices(sb)
significant_indices(sf)
