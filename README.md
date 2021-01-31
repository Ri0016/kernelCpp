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
