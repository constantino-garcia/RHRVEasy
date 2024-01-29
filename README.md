

# RHRVEasy
An R package created to automate all steps of a HRV analysis, including data preprocessing, indices calculation, and statistical analysis. The methods of this package are described in:

> García, C.A., Bardají, S., Pérez-Tirador, P., Otero, A. **RHRVEasy: heart rate variability made easy**. *Under review*

## Installation

Use `devtools` to install the package. In an R console, execute the following commands:

```R
# install.packages("devtools") # only if needed
devtools::install_github("constantino-garcia/nonlinearTseries", ref = "rqa")
devtools::install_github("constantino-garcia/RHRVEasy")
```

## Data
The folder `RRData` contains a zip file with the data used to test the package in the paper. After unzipping, these data can also be used to test the package or follow the tutorial (see next section).

## Tutorial
The `RHRVEasyTutorial.Rmd` provides a step-by-step introduction to the package. Furthermore, the results of the paper can be reproduced by completing the tutorial. To follow the tutorial, please refer to the `Data` section.


