

RHRVEasyPValues<-function(folders, correction = FALSE, correctionMethod = "bonferroni", verbose=FALSE, 
                   format = "RR", typeAnalysis = 'fourier', significance_level = 0.05, nonLinear=FALSE, ...) {
  
  dataFrameMWavelet = data.frame()
  dataFrameMTime = data.frame()
  dataFrameMFreq = data.frame()
  dataFrameMNonLinear = data.frame()
  listNonLinearStatisticalAnalysis = list()
  listTimeStatysticalAnalysis = list()
  listFreqStatysticalAnalysis = list()
  
  #We create a global variable signif_level with the level significance
  signif_level <<- significance_level
  #We create a global variable verb with the verbose mode  
  verb <<- verbose
  
  files = list()
  
  f  <- vector(mode = "list", length = length(folders))
  i=1
  
  for (folder in folders){
    file_validation(folder)
    if(i==1){
      f[[i]] =  sample(list.files(folder), size =18, replace = FALSE)
    }
    else{
      f[[i]] =  sample(list.files(folder), size =10, replace = FALSE)
    }
    dataFrameMTime = rbind(dataFrameMTime, dataFrameMTime = time_analysis(format,
                                                                          f[[i]], split_path(folder)[1], folder, ...))
    if(nonLinear == TRUE){
      cat("Performing non linear analysis...")
      dataFrameMNonLinear = rbind(dataFrameMNonLinear,non_linear_analysis(format, 
                                                                          f[[i]], split_path(folder)[1], folder, ...))
    }
    i = i+1
  }
  
  numberOfExperimentalGroups = length(folders)
  # Statistical analysis of both
  
  listTimeStatysticalAnalysis = statistical_analysisTime(dataFrameMTime,
                                                         numberOfExperimentalGroups, correctionMethod, signif_level)
  
  # FREQUENCY:
  if(typeAnalysis == "fourier"){
    i=1
    for (folder in folders){
      dataFrameMFreq = rbind(dataFrameMFreq, dataFrameMFreq = freq_analysis(format,
                                                                            f[[i]], split_path(folder)[1], folder, ...))
      i=i+1
    }
    
    listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMFreq,
                                                           numberOfExperimentalGroups, correctionMethod, signif_level)
  }
  
  # WAVELET
  if(typeAnalysis == "wavelet"){
    i=1
    for (folder in folders){
      dataFrameMWavelet = rbind(dataFrameMWavelet, dataFrameMWavelet = wavelet_analysis(format,
                                                                                        f[[i]], split_path(folder)[1], folder,
                                                                                        type = typeAnalysis, ...))
      i=i+1
    }
    
    listFreqStatysticalAnalysis = statistical_analysisFreq(dataFrameMWavelet,
                                                           numberOfExperimentalGroups, correctionMethod, signif_level)
    
    dataFrameMFreq = dataFrameMWavelet
  }
  
  
  if(!all(is.na(dataFrameMNonLinear))){
    #cat("DF non linear NOT empty: correct p values with those values")
    listNonLinearStatisticalAnalysis = statistical_analysisNonLinear(dataFrameMNonLinear,
                                                                     numberOfExperimentalGroups, correctionMethod, signif_level)
  }else{
    #cat("DF non linear empty")
    listNonLinearStatisticalAnalysis = NA
  }
  
  
  uncorrectedPvalues = colectpValues(listTimeStatysticalAnalysis, listFreqStatysticalAnalysis,
                                     listNonLinearStatisticalAnalysis)
  listpValues = correctpValues(uncorrectedPvalues, correction, correctionMethod)
  
  results = list("TimeAnalysis" = dataFrameMTime, 
                 "StatysticalAnalysisTime" = listTimeStatysticalAnalysis,
                 "FrequencyAnalysis" = dataFrameMFreq, 
                 "StatysticalAnalysisFrequency" = listFreqStatysticalAnalysis,
                 "NonLinearAnalysis" = dataFrameMNonLinear,
                 "StatysticalAnalysisNonLinear" = listNonLinearStatisticalAnalysis,
                 "pValues" = listpValues, "uncorrectedPvalues" = uncorrectedPvalues)
  
  class(results) = "RHRVEasyResult"
  results
}

ar=array(0, dim=c(10,6,14))


for(r in 1:10){
c1 = list(ULF = 0, VLF = 0, LF = 0, HF = 0,
          SDNN = 0, SDANN = 0, SDNNIDX = 0, pNN50 = 0, SDSD = 0, 
          rMSSD = 0, IRRR = 0, MADRR = 0, TINN = 0, HRVi = 0)


c2 = c1
c3 = c1
c4 = c1
c5 = c1
c6 = c1
significativityCount = list("none"=c1, "fdr"=c2, "holm"=c3, "hochberg"=c4, "BY"=c5, "bonferroni"=c6)


vec = c(significativityCount,significativityCount);

columNames = c('ULF', 'SDNN', 'IRRR', 'SDANN', 'TINN', 'HRVi', 'VLF', 'SDNNIDX',
               'LF', 'SDSD', 'rMSSD', 'MADRR', 'HF', 'pNN50')

records =c("C:\\rrsc\\normal","C:\\rrsc\\chf")

corrections = c("none", "fdr", "holm", "hochberg", "BY", "bonferroni")


  for(correct in corrections){
    for(i in 1:100){
      a=RHRVEasyPValues(folders = records,correction = TRUE, correctionMethod = correct )
      for (column in columNames){ 
        if(a$pValues[[column]] < 0.05){
           significativityCount[[correct]][[column]] = significativityCount[[correct]][[column]] +1;
        }
      }
    }
  }
  for(i in 1:6){
    for(j in 1:14){
      ar[r,i,j] = significativityCount[[i]][[j]]
    } 
  }
}

apply(ar, MARGIN = c(2,3), FUN=mean)
apply(ar, MARGIN = c(2,3), FUN=sd)

res =cbind(significativityCount$none, significativityCount$fdr,significativityCount$holm, 
           significativityCount$hochberg,significativityCount$BY,significativityCount$bonferroni)

res

