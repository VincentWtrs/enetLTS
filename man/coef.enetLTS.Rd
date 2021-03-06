\name{coef.enetLTS}
\alias{coef.enetLTS}
%-------------------------------------------------
\title{
coefficients from the \code{enetLTS} object
}
%-------------------------------------------------
\description{
  A numeric vector which extracts model coefficients from object
  returned by regression model.
}
%-------------------------------------------------
\usage{
  \method{coef}{enetLTS}(object,vers,zeros,...)
}
%-------------------------------------------------
\arguments{
\item{object}{fitted \code{enetLTS} model object.}

\item{vers}{a character string specifying for which fit to make
  predictions. Possible values are \code{reweighted} (the default) for
  predicting values from the reweighted fit, \code{raw} for predicting
  values from the raw fit.}

\item{zeros}{a logical indicating whether to keep zero coefficients
(\code{TRUE}, the default) or to omit them (\code{FALSE}).}

\item{\dots}{additional arguments from the \code{enetLTS} object if needed.}
}
%-------------------------------------------------
\value{
  a numeric vector containing the requested coefficients.
}
%-------------------------------------------------
\seealso{
\code{\link{enetLTS}},
\code{\link{predict.enetLTS}},
\code{\link{nonzeroCoef.enetLTS}}
}
%-------------------------------------------------
\examples{
## for gaussian

set.seed(86)
n <- 100; p <- 25                             # number of observations and variables
beta <- rep(0,p); beta[1:6] <- 1              # 10\% nonzero coefficients
sigma <- 0.5                                  # controls signal-to-noise ratio
x <- matrix(rnorm(n*p, sigma),nrow=n)
e <- rnorm(n,0,1)                             # error terms
eps <- 0.1                                    # contamination level
m <- ceiling(eps*n)                           # observations to be contaminated
eout <- e; eout[1:m] <- eout[1:m] + 10        # vertical outliers
yout <- c(x \%*\% beta + sigma * eout)          # response
xout <- x; xout[1:m,] <- xout[1:m,] + 10      # bad leverage points

\donttest{
fit1 <- enetLTS(xout,yout,alphas=0.5,lambdas=0.05,plot=FALSE)
coef(fit1)
coef(fit1,vers="raw")
coef(fit1,vers="reweighted",zeros=FALSE)
}


## for binomial

eps <-0.05                                     # \%10 contamination to only class 0
m <- ceiling(eps*n)
y <- sample(0:1,n,replace=TRUE)
xout <- x
xout[y==0,][1:m,] <- xout[1:m,] + 10;          # class 0
yout <- y                                      # wrong classification for vertical outliers
\dontshow{
set.seed(86)
n <- 5; p <- 15                            
beta <- rep(0,p); beta[1:6] <- 1              
sigma <- 0.5                        
x <- matrix(rnorm(n*p, sigma),nrow=n)   
e <- rnorm(n,0,1)                             # error terms
eps <- 0.1                                    # contamination level
m <- ceiling(eps*n)                           # observations to be contaminated
eout <- e; eout[1:m] <- eout[1:m] + 10        # vertical outliers
yout <- c(x \%*\% beta + sigma * eout)          # response
xout <- x; xout[1:m,] <- xout[1:m,] + 10      # bad leverage points
fit2 <- enetLTS(xout,yout,alphas=0.5,lambdas=0.05,plot=FALSE)
coef(fit2)
}
\donttest{
fit2 <- enetLTS(xout,yout,family="binomial",alphas=0.5,lambdas=0.05,plot=FALSE)
coef(fit2)
coef(fit2,vers="reweighted")
coef(fit2,vers="raw",zeros=FALSE)}
}
%-------------------------------------------------
\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
\cr Maintainer: Fatma Sevinc KURNAZ <fatmasevinckurnaz@gmail.com>;<fskurnaz@yildiz.edu.tr>}
%-------------------------------------------------
\keyword{regression}
\keyword{classification}

