

# RHRVEasy

[![Github actions CI](https://github.com/constantino-garcia/RHRVEasy/actions/workflows/github-actions.yml/badge.svg)](https://github.com/constantino-garcia/RHRVEasy/actions/workflows/github-actions.yml/badge.svg)

An R package created to automate all steps of a HRV analysis, including data preprocessing, indices calculation, and statistical analysis. The methods of this package are described in:

> García, C.A., Bardají, S., Pérez-Tirador, P., Otero, A. **RHRVEasy: heart rate variability made easy**. *Under review*

## Installation
### Installing R
1. Go to the [R project website](https://www.r-project.org/) and download the latest version of R for your operating system. In Linux systems, it may be easier to use the package manager to install R (In that case, step 2 is not necessary).
2. Install R by following the instructions provided in the website. Default options are fine for most users.

### Installing RHRVEasy
There are several options to install the package:

#### Option 1: Install from GitHub
Use `devtools` to install the package. In an R console, execute the following commands:
```R
# install.packages("devtools") # only if needed
devtools::install_github("constantino-garcia/RHRVEasy")
```

#### Option 2: Install from source 
Using the RHRVEasy_XXX.tar.gz file, where XXX is the version of the package (See [releases](https://github.com/constantino-garcia/RHRVEasy/releases/) to download it). In an R console, execute the following commands:
```R
install.packages("path/to/RHRVEasy_XXX.tar.gz", repos = NULL, type = "source")
```

#### Troubleshooting
In case dependencies are not installed automatically, you can install them manually by running:
```R
install.packages(c("boot", "broom", "doSNOW", "foreach",
  "iterators", "nonlinearTseries", "parallel", "plotrix",
  "PMCMRplus", "progress", "RHRV", "segmented", 
  "tibble", "tidyr", "writexl"))
```


## API overview 

The main function of the package is `RHRVEasy` and takes a single mandatory argument:
a list of folders, each containing the recordings of a same 
population. 

```R 
easyAnalysis <- RHRVEasy(c("path/to/folder1", "path/to/folder2")) 
```

`RHRVEasy` calculates time, frequency, and nonlinear domain HRV indices, 
and then it applies hypothesis test, and corrects the significance levels. If 
there are more than two experimental groups and statistically significant 
differences are found, it performs a post-hoc analysis to find out which groups 
have the differences. 

More details about the API can be found in the package documentation `RHRVEasy_XXX.pdf` and the tutorial `RHRVEasyTutorial.R` (or its compiled version `RHRVEasyTutorial.pdf`).


## Data
The folder `RRData` contains a zip file with the data used to test the package in the paper. After unzipping, these data can also be used to test the package or follow the tutorial (see next section).

## Tutorial
The `RHRVEasyTutorial.Rmd` provides a step-by-step introduction to the package. Furthermore, the results of the paper can be reproduced by completing the tutorial. To follow the tutorial, please refer to the `Data` section and unzip the zip folder under the `RRData` directory.


