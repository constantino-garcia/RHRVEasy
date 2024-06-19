fileValidation <- function(path, easyOptions) {
  # 1. Check if path really exists
  if (dir.exists(path) != TRUE) {
    stop("The path ", path, " does not exist")
  }
  # 2. The path contains files:
  if ((length(list.files(path)) > 0) != TRUE) {
    stop(
      paste0(
        "The path ", path, " exists but there are no files in it"
      )
    )
  } else {
    if (easyOptions$verbose) {
      message(
        "The path ", path, " exists and there are files in it"
      )
    }
  }
}


#' @importFrom RHRV CreateHRVData SetVerbose LoadBeat BuildNIHR FilterNIHR
prepareAnalysis <- function(file, rrs, format, easyOptions) {
  hrv.data <- CreateHRVData()
  hrv.data <- SetVerbose(hrv.data, FALSE)

  hrv.data <- tryCatch({
    hrv.data <- LoadBeat(
      fileType = format,
      HRVData = hrv.data,
      Recordname = file,
      RecordPath = rrs
    )
    hrv.data
  },
  error = function(cond) {
    stop(
      paste(
        "The file \"",
        file,
        "\" could not be loaded. Check if the file is in the correct format; the specified format was \"",
        format,
        "\".",
        sep = ""
      )
    )
  })
  hrv.data <- withCallingHandlers(
    {
      hrv.data <- BuildNIHR(hrv.data)
      hrv.data <- FilterNIHR(hrv.data)
      hrv.data$Beat <- hrv.data$Beat[2:nrow(hrv.data$Beat), ]
      hrv.data
    },
    warning = function(w) {
      w$message <- paste0("File ", file,": ", w$message, "\n")
      warning(w)
      invokeRestart("muffleWarning")
    }
  )
  hrv.data
}

#Calls an RHRV function with hrv.data after cleaning the parameters
#' @importFrom plotrix clean.args
easyCall <- function(hrv.data, mf, ...) {
  args.list <- plotrix::clean.args(list(...), mf)
  args.list$HRVData <- hrv.data
  do.call(mf, args.list)
}

#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom doSNOW registerDoSNOW
prepareEasyCluster <- function(nJobs, verbose, clusterLogFile) {
  nCores <- parallel::detectCores(logical = FALSE)
  if (nJobs <= 0) {
    nJobs <- nCores
  } else if (nJobs > nCores) {
    nJobs <- nCores
  }
  if (nJobs > 1) {
    cl <- parallel::makeCluster(nJobs, outfile=clusterLogFile) # using outfile = "" may be useful for debugging
    doParallel::registerDoParallel(cl)
    #doSNOW::registerDoSNOW(cl)
    if (verbose) {
      message(paste("Registering cluster with", nJobs, "nodes"))
    }
  } else {
    cl <- NULL
    foreach::registerDoSEQ()
  }
  cl
}


splitPath <- function(path) {
  if (dirname(path) %in% c(".", path))
    return(basename(path))
  return(c(basename(path), splitPath(dirname(path))))
}




# FIXME @importClassesFrom progress progress_bar

#' @import progress
updateProgressFactory <- function(analysis, files){
  pb <- progress::progress_bar$new(
    format = paste(analysis, "of :file [:bar] :elapsed | eta: :eta"),
    total = length(files),
    width = 80
  )
  function(n) {
    pb$tick(tokens = list("file" = files[n]))
  }
}
