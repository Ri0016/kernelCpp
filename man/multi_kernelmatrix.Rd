\name{multi_kernelmatrix}
\alias{multi_kernelmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{multi_kernelmatrix
%%  ~~function to do ... ~~
}
\description{
return a matrix of kernel density of input matrix X.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
multi_kernelmatrix(X,h,kernel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{input matrix
%%     ~~Describe \code{x} here~~
} \item{h}{bandwidth
} \item{kernel}{kernel function type

Triangle:kernel = "triangle" or "tria";

Epanechnikov:kernel = "epanechnikov" or "e";

Quartic:kernel = "quartic" or "q";

Triweight:kernel = "triweight" or "triw";

Gaussian:kernel = "gaussian" or "g";

Cosine:kernel = "cosine" or "c";
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
If you want to remove the diag of result matrix, following the code:

diag(result) <- 0

where the result is the function return value, like:

result = multi_kernelmatrix(X,h,"c")

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The return value is a symmetric matrix. The (i,j) element is the product of kernel density  in each dimension between row i and row j of input matrix.

}
\references{
%% ~put references to the literature/web site here ~

Univariate kernel demo:

https://mathisonian.github.io/kde/


Multivariate kernel density estimation:

https://bookdown.org/egarpor/NP-UC3M/kde-ii-mult.html


Different kernel functions:

https://en.wikipedia.org/wiki/Kernel_(statistics)

}
\author{
%%  ~~who you are~~
Ri Wang
}
\note{
%%  ~~further notes~~
remove the package with following command:

remove.packages("kernelCpp",.libPaths())

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
function (x)
{
  }
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ kerneldensity }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
