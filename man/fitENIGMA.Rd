\name{fitENIGMA}
\alias{fitENIGMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of unigram models for metagenome analysis
}
\description{
Estimate the parameteres regressive unigram mixtutre model using Laplace approximation.
}
\usage{
fitENIGMA(Y, formula = NULL, X = NULL, L = 2, alpha = NULL, data = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Matrix of multinomial response.
}
  \item{formula}{
Formula of linear predictor.
}
  \item{X}{
Design matrix (instead of formula).
}
  \item{L}{
Number of components (Numeber of cluster).
}
  \item{alpha}{
Hyperparameters of dirichlet prior (default 1)
}
  \item{data}{
An optional data frame which containing the variables in the model.
}
  \item{\dots}{
additional arguments to be passed to optimizing function in rstan.
}
}
\value{

\item{pi }{mixting ratio}
\item{beta }{regression coefficients}
\item{gamma }{The baseline parameter (intercept) which allowed change by clusters}
\item{sigma }{The hyperparameter of beta}
\item{tau }{The hyper parameter of gamma}
\item{cate }{The probability which a sample belong each class}
\item{Xbeta }{The linear predictor}
}
\references{
Ko Abe, Masaaki Hirayama, Kinji Ohno, and Teppei Shimamura, An Enterotype-like Unigram Class Model for Identifying Microbial Associations with Diseases, submitted.
}
\author{
Ko ABE
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
softmax <-function(x){
  tmp = max(x)
  return(exp(x-tmp)/sum(exp(x-tmp)))
}

# generate synthetic data
N <- 100
L <- 3
K <- 100
set.seed(1234)
gamma <- matrix(rnorm(L*K),L,K)
beta <- rnorm(K)
X <- matrix(rbinom(N,1,0.5))
z <-sample.int(L,N,replace = TRUE)
p <-t(apply(gamma[z,] + (X \%*\% beta),1,softmax))
y <-t(apply(p,1,function(x)rmultinom(1,2000,x)))

fit <- fitENIGMA(y,X=X,L = 3)
plot(beta,fit$par$beta)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
