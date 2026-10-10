# Run from an upstream STA checkout after applying issue-6.patch:
# Rscript STACMR-R/tests/test-mr-bootstrap-pars.R
# Does not install packages, edit the checkout, or run staCMRsetup.R.
library(rJava)
library(methods)
source("STACMR-R/Utility functions/gen2list.R")
source("STACMR-R/java/jMRfits.R")
source("STACMR-R/staMRFIT.R")
.jinit(classpath = normalizePath("STACMR-R/java/fxMR-0.3.44.jar"))
delay <- read.table("STACMR-R/Data files/delay.dat")
order <- list(1:4, 5:8, c(5, 1), c(6, 2), c(7, 3), c(8, 4))
out <- staMRFIT(delay, partial = order, nsample = 4)
stopifnot(length(out$pars) == 2L, length(out$fits) == 4L,
          abs(out$datafit - 0.1721) < 5e-5)
for (p in out$pars) stopifnot(identical(dim(p), c(4L, 8L)), all(is.finite(p)))
cat("Issue #6 regression passed: both DVs have 4 x 8 bootstrap estimates.\n")
