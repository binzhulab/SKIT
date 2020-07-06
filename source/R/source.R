skit <- function(x, y, bandwidth=NULL, nboot=1000, print=0) {

  if (!is.vector(x)) stop("ERROR: x must be a vector")
  if (!is.vector(y)) stop("ERROR: y must be a vector")
  if (length(x) != length(y)) stop("ERROR: x and y should have the same length")
  tmp <- is.finite(x) & is.finite(y)
  if (!all(tmp)) {
    x <- x[tmp]
    y <- y[tmp]
    if (!length(x)) stop("ERROR: all observations have been removed")
  }
  if (nboot < 0) nboot <- 0
  bwFlag <- length(bandwidth)
  if (bwFlag) {
    if ((bwFlag > 1) || (bandwidth < 0) || !is.finite(bandwidth)) stop("ERROR: invalid value for bandwidth")
    sigma    <- NULL
  } 

  ret <- ts_main(x, y, nboot, bandwidth, bwFlag, print)

  ret

} # END: skit

getDefaultBandwidth <- function(x, y, default=1) {

  n    <- length(x)
  x1   <- x[x != 0]
  y1   <- y[y != 0]
  sigx <- sd(x1)
  sigy <- sd(y1)
  SIG  <- ifelse(is.na(sigx) & is.na(sigy), 1, max(c(sigx,sigy), na.rm = TRUE))
  if(SIG == 0) bw <- default else bw <- SIG*n^(-0.2)

  list(bandwidth=bw, sigma=SIG)

} # END: getDefaultBandwidth

ts_main <- function(x, y, nboot, bnd, bwFlag, print) {

  n     <- length(x)
  tests <- rep(-9999, 5)
  
  if (!nboot) {
    if (!bwFlag) bnd <- getDefaultBandwidth(x, y)$bandwidth 
    tmp   <- .C("C_TS", as.numeric(x), as.numeric(y), as.integer(n), as.numeric(bnd),
                 ret=as.numeric(tests), PACKAGE="SKIT")
    tests <- tmp$ret
    pvals <- rep(NA, length(tests))
  } else {  
    pvals  <- rep(-9999, 5)
    retbnd <- -9999
    if (!bwFlag) bnd <- -9999
    tmp   <- .C("C_TS_boot", as.numeric(x), as.numeric(y), as.integer(n), 
               as.numeric(bnd), as.integer(bwFlag), as.integer(nboot), 
               as.integer(print), retp=as.numeric(pvals), retbw=as.numeric(retbnd), 
               retobs=as.numeric(tests),
               PACKAGE="SKIT")
    tests <- tmp$retobs
    pvals <- tmp$retp
    bnd   <- tmp$retbw
  } 
  names(tests)  <- names(pvals) <- c("T", "T1", "T2", "T3", "T4")

  list(tests=tests, pvalues=pvals, bandwidth.obs=bnd)

} # END: ts_main


