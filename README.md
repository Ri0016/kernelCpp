# kernelCpp

利用RcppParallel,计算核密度矩阵.

## 环境要求:

需要电脑中有c++环境, 在R中安装 Rcpp, RcppParallel包.

install.packages("Rcpp")

install.packages("RcppParallel")

## 安装:

1.  利用devtools包,自动下载安装

install.packages("devtools")

在R控制台输入命令: devtools::install_github("Ri0016/kernelCpp",force =
TRUE)

2.  手动下载zip,通过文件离线安装

## 使用

library("kernelCpp")

查看函数说明:

?multi_kernelmatrix



## 卸载

remove.packages("kernelCpp",.libPaths())

## 相关链接:

RcppParallel: <https://rcppcore.github.io/RcppParallel/>

一维核估计动态演示: <https://mathisonian.github.io/kde/>

多维核估计: <https://bookdown.org/egarpor/NP-UC3M/kde-ii-mult.html>

不同核函数类别: <https://en.wikipedia.org/wiki/Kernel_(statistics)>



## 验证

```r
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
p=20
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
```
