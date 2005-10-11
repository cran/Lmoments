\name{Lmoments}
\alias{Lmoments}
\alias{Lcoefs}
\alias{Lmomcov}
\concept{moments}
\concept{skewness}
\concept{kurtosis}
\concept{descriptive statistics}
\title{L-moments}
\description{
  Calculates sample L-moments, L-coefficients and covariance matrix of L-moments. 
}
\usage{
Lmoments(data,rmax=4)
Lcoefs(data)
Lmomcov(data,rmax=4)
}
\arguments{
  \item{data}{matrix or data frame}
  \item{rmax}{maximum order of L-moments}
}
\value{\code{Lmoments} returns an array of L-moments containing a row for each variable in data. 
\code{Lcoefs} returns an array of L-coefficients (L-mean, L-scale, L-skewness, L-kurtosis) 
containing a row for each variable in data.
Lmomcov returns the covariance matrix of L-moments.}
\references{Karvanen, J. 2005. Estimation of quantile mixtures via L-moments and trimmed L-moments, 
Accepted for publication in \emph{Computational Statistics & Data Analysis} 
\url{http://www.bsp.brain.riken.jp/publications/2005/karvanen_quantile_mixtures.pdf}.

Elamir, E. A., Seheult, A. H. 2004. Exact variance structure of sample L-moments, \emph{Journal of 
Statistical Planning and Inference} \bold{124} (2) 337--359.

Hosking, J. 1990. L-moments: Analysis and estimation distributions using linear combinations of order statistics, 
\emph{Journal of Royal Statistical Society B} \bold{52}, 105--124.
}
\seealso{\code{\link{t1lmoments}} for trimmed L-moments, and 
\code{\link{dnormpoly}}, \code{\link{lmom2normpoly4}} and \code{\link{covnormpoly4}} for the normal-polynomial quantile mixture
}
\examples{
#Generates a sample 500 observations from the normal-polynomial quantile mixture, 
#calculates the L-moments and their covariance matrix,
#estimates parameters via L-moments and 
#plots the true pdf and the estimated pdf together with the histogram of the data.
true_params<-lmom2normpoly4(c(0,1,0.2,0.05));
x<-rnormpoly(500,true_params);
lmoments<-Lmoments(x);
lmomcov<-Lmomcov(x);
estim_params<-lmom2normpoly4(lmoments);
hist(x,30,freq=FALSE)
plotpoints<-seq(min(x)-1,max(x)+1,by=0.01);
lines(plotpoints,dnormpoly(plotpoints,estim_params),col='red');
lines(plotpoints,dnormpoly(plotpoints,true_params),col='blue');
}
\author{Juha Karvanen <\email{juha.karvanen@ktl.fi}>}
\keyword{robust}
\keyword{univar}

