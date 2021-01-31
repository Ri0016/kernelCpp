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

\examples{
##compare with functions implemented in R

Epanechnikov <- function (x)
{

  xPh<- abs(x)
  xPh[xPh <=1] <-1
  xPh[xPh>1] <- 0
  kernalX <- 0.75*(1-x^2)*xPh
  return(kernalX)
}

k_kernel<-function(h,x){
  #h smoothing parameters
  #X regressors
  p=dim(x)[2]
  temp=1
  for(j in 1:p){
    xx=outer(x[,j],x[,j],'-')
    k_value=Epanechnikov(xx/h)
    temp=temp*k_value/h
  }
  temp=temp#-diag(diag(temp))
  return(temp)
}
epan_kernel<-function(h,x){
  #h smoothing parameters
  #X regressors
  p=dim(x)[2]
  temp=1
  for(j in 1:p){
    xx=outer(x[,j],x[,j],'-')
    k_value=0.75*(1-(xx/h)^2)*(abs(xx/h) < 1)
    temp=temp*k_value/h
  }
  temp=temp#-diag(diag(temp))
  return(temp)
}
quartic_kernel<-function(h,x){
  #h smoothing parameters
  #X regressors
  p=dim(x)[2]
  temp=1
  for(j in 1:p){
    xx=outer(x[,j],x[,j],'-')
    k_value=0.9375*(1-(xx/h)^2)^2*(abs(xx/h) < 1)
    temp=temp*k_value/h
  }
  temp=temp#-diag(diag(temp))
  return(temp)
}
triweight_kernel<-function(h,x){
  #h smoothing parameters
  #X regressors
  p=dim(x)[2]
  temp=1
  for(j in 1:p){
    xx=outer(x[,j],x[,j],'-')
    k_value=1.09375*(1-(xx/h)^2)^3*(abs(xx/h) < 1)
    temp=temp*k_value/h
  }
  temp=temp#-diag(diag(temp))
  return(temp)
}
triangle_kernel<-function(h,x){
  #h smoothing parameters
  #X regressors
  p=dim(x)[2]
  temp=1
  for(j in 1:p){
    xx=outer(x[,j],x[,j],'-')
    k_value=(1-abs(xx/h))*(abs(xx/h) < 1)
    temp=temp*k_value/h
  }
  temp=temp#-diag(diag(temp))
  return(temp)
}

cosine_kernel<-function(h,x){
  #h smoothing parameters
  #X regressors
  p=dim(x)[2]
  temp=1
  for(j in 1:p){
    xx=outer(x[,j],x[,j],'-')
    k_value=pi/4*cos(pi*(xx/h)/2)*(abs(xx/h) < 1)
    temp=temp*k_value/h
  }
  temp=temp#-diag(diag(temp))
  return(temp)
}

normal_kernel<-function(h,x){
  #h smoothing parameters
  #X regressors
  p=dim(x)[2]
  temp=1
  for(j in 1:p){
    xx=outer(x[,j],x[,j],'-')
    k_value=dnorm(xx/h)
    temp=temp*k_value/h
  }
  temp=temp#-diag(diag(temp))
  return(temp)
}
n=300
p=5
h =0.1
X <- matrix(runif(n*p),n,p)
e <- rnorm(n)
#epan
k_kernel(h,X) -> r1
multi_kernelmatrix(X,h,"e") -> r2
epan_kernel(h,X) -> r3
bench::mark(
  k_kernel(h,X),
  epan_kernel(h,X),
  multi_kernelmatrix(X,h,"e"),
  check = TRUE,
  relative = TRUE
)

all.equal(r1,r2)
all.equal(r1,r3)

#norm
normal_kernel(h,X) -> t1
multi_kernelmatrix(X,h,"g") -> t2
bench::mark(
  normal_kernel(h,X),
  multi_kernelmatrix(X,h,"g"),
  check = TRUE,
  relative = TRUE
)
all.equal(t1,t2)

#quartifc

quartic_kernel(h,X) ->q1
multi_kernelmatrix(X,h,"q") -> q2
bench::mark(
  quartic_kernel(h,X),
  multi_kernelmatrix(X,h,"q"),
  check = TRUE,
  relative = TRUE
)
all.equal(q1,q2)

#triweight

triweight_kernel(h,X) -> tw1
multi_kernelmatrix(X,h,"triw") -> tw2
bench::mark(
  triweight_kernel(h,X),
  multi_kernelmatrix(X,h,"triw") ,
  check = TRUE,
  relative = TRUE
)
all.equal(tw1,tw2)

#triangle
triangle_kernel(h,X) -> ta1
multi_kernelmatrix(X,h,"tria") -> ta2
bench::mark(
  triangle_kernel(h,X),
  multi_kernelmatrix(X,h,"tria") ,
  check = TRUE,
  relative = TRUE
)
all.equal(ta1,ta2)

#cosine_kernel

cosine_kernel(h,X) -> c1
multi_kernelmatrix(X,h,"c") -> c2
bench::mark(
  cosine_kernel(h,X),
  multi_kernelmatrix(X,h,"c"),
  check = TRUE,
  relative = TRUE
)
all.equal(c1,c2)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ kerneldensity }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
