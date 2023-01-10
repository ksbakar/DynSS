
################################################################################
################################################################################

## implementation file for DynSS method
## to reproduce the results in the paper

source("source_file_DynSS.R")


fig_1()
## Figure 2
fig_2()
## Figure 3
fig_3()
## Figure 4
fig_4()
## Figure 5
fig_5()
## Figure 6
fig_6()
## Figure 7
fig_7(margin_of_error = 0.1) # 
fig_7(margin_of_error = 0.05) # not reported in the paper
##


################################################################################
################################################################################

## Code to reproduce results in the Supplementary section

## Figure S2
fig_S2(D=list(c(0.05,0.95)))
fig_S2(D=list(c(0.01,0.99)))
fig_S2(D=list(c(0.1,0.9))) # not included in the paper

## Figure S3
fig_S3(D=list(c(0.05,0.95)),margin_of_error=0.05)
fig_S3(D=list(c(0.1,0.9)))
fig_S3(D=list(c(0.01,0.99)))

## Figure S4
fig_S4 <- function(D0=0.025,D1=0.975){
  fig_6(drop_arm_type=c("futility"),D0=D0,D1=D1)
}
fig_S4(D0=0.1,D1=0.9)
fig_S4(D0=0.05,D1=0.95)
fig_S4(D0=0.01,D1=0.99)

## Figure S5
fig_S5 <- function(D0=0.025,D1=0.975,margin_of_error=0.10){
  fig_7(drop_arm_type=c("futility"),D0=D0,D1=D1,margin_of_error = margin_of_error)
}
fig_S5(D0=0.1,D1=0.9)
fig_S5(D0=0.05,D1=0.95)
fig_S5(D0=0.01,D1=0.99)

## Figure S6
fig_S6 <- function(D0=0.025,D1=0.975){
  fig_6(drop_arm_type=c("efficacy"),D0=D0,D1=D1)
}
fig_S6(D0=0.1,D1=0.9)
fig_S6(D0=0.05,D1=0.95)
fig_S6(D0=0.01,D1=0.99)

## Figure S7
fig_S7 <- function(D0=0.025,D1=0.975,margin_of_error=0.10){
  fig_7(drop_arm_type=c("efficacy"),D0=D0,D1=D1,margin_of_error = margin_of_error)
}
fig_S7(D0=0.1,D1=0.9)
fig_S7(D0=0.05,D1=0.95)
fig_S7(D0=0.01,D1=0.99)


################################################################################
################################################################################