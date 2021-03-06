\name{plotResid.enetLTS}
\alias{plotResid.enetLTS}
%-------------------------------------------------
\title{
residuals plots from the \code{"enetLTS"} object
}
%-------------------------------------------------
\description{
  Produce plots for the residuals of the current model.
}
%-------------------------------------------------
\usage{
plotResid.enetLTS(object,vers=c("reweighted","raw"), ...)
}
%------------------------------------
\arguments{
\item{object}{the model fit to be plotted.}

\item{vers}{a character string denoting which model to use for the plots.
  Possible values are \code{"reweighted"} (the default) for
  plots from the reweighted fit, and \code{"raw"} for
  plots from the raw fit.}

\item{\dots}{additional arguments from the \code{enetLTS} object if needed.}
}
%-------------------------------------------------
\value{
  An object of class \code{"ggplot"} (see \code{\link[ggplot2]{ggplot}}).
}
%-------------------------------------------------
\note{
gives the plot of
        - residuals vs indices. (for both \code{family="binomial"} and \code{family="gaussian"}).

        - additionally, residuals vs fitted values (for only \code{family="gaussian"}).
}
%-------------------------------------------------
\seealso{
\code{\link[ggplot2]{ggplot}},
\code{\link{enetLTS}},
\code{\link{predict.enetLTS}},
\code{\link{residuals.enetLTS}},
\code{\link{fitted.enetLTS}}
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
plotResid.enetLTS(fit1)
plotResid.enetLTS(fit1,vers="raw")}

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
plotResid.enetLTS(fit2)
}
\donttest{
fit2 <- enetLTS(xout,yout,family="binomial",alphas=0.5,lambdas=0.05,plot=FALSE)
plotResid.enetLTS(fit2)
plotResid.enetLTS(fit2,vers="raw")}
}
%-------------------------------------------------
\references{
Kurnaz, F.S., Hoffmann, I. and Filzmoser, P. (2017) Robust and sparse
estimation methods for high dimensional linear and logistic regression.
\emph{Chemometrics and Intelligent Laboratory Systems}.
}
%-------------------------------------------------
\author{
Fatma Sevinc KURNAZ, Irene HOFFMANN, Peter FILZMOSER
\cr Maintainer: Fatma Sevinc KURNAZ <fatmasevincskurnaz@gmail.com>;<fskurnaz@yildiz.edu.tr>}
%-------------------------------------------------
\keyword{regression}
\keyword{classification}


