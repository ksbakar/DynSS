
################################################################################
################################################################################

## implementation file for DynSS method
## to reproduce the results in the paper

library(devtools)
library(roxygen2)
## Needed because this file has roxygen2 comments. Otherwise you get a
## 'could not find function 'digest'' error
source_url("https://github.com/ksbakar/DynSS/blob/main/source_file_DynSS.R?raw=TRUE")
source("source_file_DynSS.R")

package_list <- c("dplyr", "data.table", "spTimer", "LaplacesDemon", "extraDistr",
                  "pbapply", "MASS", "ggplot2", "grid", "gridExtra", "ggridges")
checkPkg <- function(package_list){
  check_installed <- package_list[!(package_list %in% installed.packages()[ , "Package"])]
  if(length(check_installed)) {
    print("Installing required packages from CRAN ...")
    install.packages(check_installed)
    print("All required packages are installed in the working machine ...")
  }
  else{
    print("All required packages are installed in the working machine ...")
  }
}
checkPkg(package_list)
lapply(package_list, require, character.only = TRUE)


## Figure 1
fig_1()
## Figure 2
fig_2()
## Figure 3
fig_3()
## Figure 4
fig_4()
## Figure 5
fig_5()
##


################################################################################
################################################################################

## Code to reproduce results in the Supplementary section

## Fig S1
fig_S1()
## Fig S3
fig_S3(D=list(c(0.05,0.95)))
fig_S3(D=list(c(0.05,0.99)))
## Fig S4
fig_S4(margin_of_error=0.05)
fig_S4(margin_of_error=0.01)
## Fig S5
fig_4(drop_arm_type=c("futility"),D0=0.1,D1=0.9)
fig_4(drop_arm_type=c("futility"),D0=0.05,D1=0.95)
fig_4(drop_arm_type=c("futility"),D0=0.01,D1=0.99)
## Fig S6
fig_5(drop_arm_type=c("futility"),D0=0.1,D1=0.9)
fig_5(drop_arm_type=c("futility"),D0=0.05,D1=0.95)
fig_5(drop_arm_type=c("futility"),D0=0.01,D1=0.99)
## Fig S7
fig_4(drop_arm_type=c("efficacy"),D0=0.1,D1=0.9)
fig_4(drop_arm_type=c("efficacy"),D0=0.05,D1=0.95)
fig_4(drop_arm_type=c("efficacy"),D0=0.01,D1=0.99)
## Fig S8
fig_5(drop_arm_type=c("efficacy"),D0=0.1,D1=0.9)
fig_5(drop_arm_type=c("efficacy"),D0=0.05,D1=0.95)
fig_5(drop_arm_type=c("efficacy"),D0=0.01,D1=0.99)

################################################################################
################################################################################
