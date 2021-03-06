\name{calcSE}
\alias{calcSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate standard error from an object made by fitENIGMA.
}
\description{
Calculate standard error from inverse Hessian matrix
}
\usage{
calcSE(fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
enter the object made by the function fitENIGMA.
}
}
\value{
A numeric vector with the parameteres standard error.
}
\references{
%% ~put references to the literature/web site here ~
Ko Abe, Masaaki Hirayama, Kinji Ohno, and Teppei Shimamura, An Enterotype-like Unigram Class Model for Identifying Microbial Associations with Diseases, submitted.
}
\author{
Ko ABE
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
calcSE(fit)
}
