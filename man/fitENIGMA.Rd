\name{fitENIGMA}
\alias{fitENIGMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fitENIGMA(Y, formula = NULL, X = NULL, L = 2, alpha = NULL, data = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y, formula = NULL, X = NULL, L = 2, alpha = NULL, data = NULL, 
    ...) 
{
    if (is.null(X)) {
        X <- model.matrix(formula, data = data)
    }
    if (is.null(alpha)) {
        alpha = rep(1, L)
    }
    dat4stan <- list(N = nrow(Y), K = ncol(Y), D = ncol(X), L = L, 
        y = Y, X = X, alpha = alpha)
    rstan::optimizing(stanmodels$softmaxreg, dat4stan, as_vector = FALSE, 
        hessian = TRUE, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
