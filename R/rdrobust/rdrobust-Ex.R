pkgname <- "rdrobust"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('rdrobust')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("rdbwselect")
### * rdbwselect

flush(stderr()); flush(stdout())

### Name: rdbwselect
### Title: Bandwidth Selection Procedures for Local Polynomial Regression
###   Discontinuity Estimators
### Aliases: rdbwselect print.rdbwselect summary.rdbwselect

### ** Examples

x<-runif(1000,-1,1)
y<-5+3*x+2*(x>=0)+rnorm(1000)
rdbwselect(y,x)



cleanEx()
nameEx("rdbwselect_2014")
### * rdbwselect_2014

flush(stderr()); flush(stdout())

### Name: rdbwselect_2014
### Title: Deprecated Bandwidth Selection Procedures for Local-Polynomial
###   Regression-Discontinuity Estimators.
### Aliases: rdbwselect_2014 print.rdbwselect_2014 summary.rdbwselect_2014

### ** Examples

x<-runif(1000,-1,1)
y<-5+3*x+2*(x>=0)+rnorm(1000)
rdbwselect_2014(y,x)



cleanEx()
nameEx("rdplot")
### * rdplot

flush(stderr()); flush(stdout())

### Name: rdplot
### Title: Data-Driven Regression Discontinuity Plots
### Aliases: rdplot print.rdplot summary.rdplot
### Keywords: regression discontinuity RD plots binning partitioning tuning
###   parameter selection

### ** Examples

x<-runif(1000,-1,1)
y<-5+3*x+2*(x>=0)+rnorm(1000)
rdplot(y,x)



cleanEx()
nameEx("rdrobust")
### * rdrobust

flush(stderr()); flush(stdout())

### Name: rdrobust
### Title: Local-Polynomial RD Estimation with Robust Confidence Intervals
### Aliases: rdrobust print.rdrobust summary.rdrobust
### Keywords: RDD Robust Estimation

### ** Examples
 
x<-runif(1000,-1,1)
y<-5+3*x+2*(x>=0)+rnorm(1000)
rdrobust(y,x)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
