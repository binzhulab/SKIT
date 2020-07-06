
\name{skit}
\alias{skit}
\title{Semiparametric Kernel Independence Test (SKIT) with excess zeros}
\description{
Performs the test of independence between two vectors when there are excess zeros.
}
\usage{ skit(x, y, bandwidth=NULL, nboot=1000, print=0)
}
\arguments{
  \item{x}{Numeric vector possibly containing zeros.}
  \item{y}{Numeric vector possibly containing zeros with same length as \code{x}.}
  \item{bandwidth}{Bandwith parameter (see details). The default is NULL so that a default
          value will be computed from the data.}
  \item{nboot}{Number of bootstrap samples to estimate p-values. If p-values
               are not desired, then set to 0.
            The default is 1000.}
  \item{print}{0 or 1 to print information. The default is 0.}
}
\details{
  Non-numeric values in \code{x} or \code{y} will be removed from both vectors.

  If \code{bandwidth = NULL}, then the \code{bandwidth} parameter will be computed 
  for the observed data, and re-computed for each bootstrap sample.
  It is computed as \code{sigma*n^(-0.2)}, where \code{n = length(x)},
  \code{sigma = max(c(sd(x[x != 0]), sd(y[y != 0])))}.

}
\value{
  A list containing the observed test statistics, bootstrap-estimated p-values,
  and the bandwith parameter for the observed data.

  The overall test statistic \eqn{\widehat{T}} is computed as \eqn{\widehat{T}_1 + \widehat{T}_2 + \widehat{T}_3 + \widehat{T}_4} 
  and for each component, p-value is estimated via bootstrap.
}
%\references{}
%\author{Bin Zhu <bin.zhu@nih.gov>, Lisa Mirabello and Nilanjan Chatterjee}
\examples{
  set.seed(123)

  n   <- 500
  b00 <- 0.09; b10 <- 0.21; b01 <- 0.21; b11 <- 0.49
  id  <- sample(1:4, n, replace = TRUE, prob = c(b00, b10, b01, b11))
  subn <- table(factor(id, levels = 1:4))
  n2 <- subn[names(subn) == "2"]
  n3 <- subn[names(subn) == "3"]
  n4 <- subn[names(subn) == "4"]
  Biv <- matrix(0, nrow = n, ncol = 2)
  if(n2 != 0) Biv[id == 2,1] <- rnorm(n2, 0.4, 0.1)
  if(n3 != 0) Biv[id == 3,2] <- rnorm(n3, 0.4, 0.1)
  if(n4 != 0) Biv[id == 4,] <- cbind(rnorm(n4, 0.4, 0.1), rnorm(n4, 0.4, 0.1))
  x <- Biv[,1]
  y <- Biv[,2]

  skit(x, y, nboot=0)

}
\keyword{ test } 
%\seealso{\code{\link{REBET}}}
