#include <Rcpp.h>
#include <RcppParallel.h>
#include <math.h>

using namespace Rcpp;
using namespace RcppParallel;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::depends(RcppParallel)]]








//RcppParallel without boost
template <typename InputIterator1, typename InputIterator2>
inline double Epan_kernel(InputIterator1 begin1,InputIterator2 begin2,
                          InputIterator1 end1, double h) {

  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  double prod=1.0;
  while(it1 != end1){
    double d1 = *it1++;
    double d2 = *it2++;
    double tx =pow((d1-d2)/h,2);
    double tp =0;
    if ( tx >= 1 ) {
      tp = 0;
    }else{
      tp =0.75*(1- tx)/h ;
    }
    prod=prod*tp;
  }

  return prod;
}
struct KernelEpanechnikov : public Worker
{
  // source matrix
  const RMatrix<double> input;

  double h;
  int p;
  // destination matrix
  RMatrix<double> output;

  // initialize with source and destination
  KernelEpanechnikov(NumericMatrix input, NumericMatrix output, double h )
    : input(input),h(h), output(output) {}


  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++){
      for (std::size_t j = 0; j <= i; j++){
        //prod = 1.0;
        RMatrix<double>::Row row1 = input.row(i);
        RMatrix<double>::Row row2 = input.row(j);
        // function to apply
        output(i,j) = Epan_kernel(row1.begin(),row2.begin(),row1.end(),h);
        output(j,i) = output(i,j);
      }

    }
  }
};

NumericMatrix KernelEpanechnikov_RcppParallel(NumericMatrix x, double h) {


  NumericMatrix rmat(x.nrow(), x.nrow());

  // declare the  instance
  KernelEpanechnikov Kh_epan(x,rmat,h);

  // call parallelFor to start the work
  parallelFor(0, x.nrow(), Kh_epan);

  // return the computed matrix
  return rmat;
}







// gaussian without boost
template <typename InputIterator1, typename InputIterator2>
inline double Gaussian_kernel(InputIterator1 begin1,InputIterator2 begin2,
                              InputIterator1 end1, double h) {

  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  double prod=1.0;
  while(it1 != end1){
    double d1 = *it1++;
    double d2 = *it2++;
    double tp = R::dnorm((d1-d2)/h,0,1,false)/h;
    prod=prod*tp;
  }

  return prod;
}
struct GaussianKernel : public Worker
{
  // source matrix
  const RMatrix<double> input;

  double h;
  int p;
  // destination matrix
  RMatrix<double> output;

  // initialize with source and destination
  GaussianKernel(NumericMatrix input, NumericMatrix output, double h )
    : input(input),h(h), output(output) {}


  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++){
      for (std::size_t j = 0; j <= i; j++){
        //prod = 1.0;
        RMatrix<double>::Row row1 = input.row(i);
        RMatrix<double>::Row row2 = input.row(j);
        // function to apply
        output(i,j) = Gaussian_kernel(row1.begin(),row2.begin(),row1.end(),h);
        output(j,i) = output(i,j);
      }

    }
  }
};

NumericMatrix KernelGaussian_RcppParallel(NumericMatrix x, double h) {


  NumericMatrix rmat(x.nrow(), x.nrow());

  // declare the  instance
  GaussianKernel Kh_epan(x,rmat,h);

  // call parallelFor to start the work
  parallelFor(0, x.nrow(), Kh_epan);

  // return the computed matrix
  return rmat;
}




// quartic without boost
template <typename InputIterator1, typename InputIterator2>
inline double Quartic_kernel(InputIterator1 begin1,InputIterator2 begin2,
                             InputIterator1 end1, double h) {

  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  double prod=1.0;
  while(it1 != end1){
    double d1 = *it1++;
    double d2 = *it2++;
    double tx =pow((d1-d2)/h,2);
    double tp =0;
    if ( tx >= 1 ) {
      tp = 0;
    }else{
      tp =0.9375*pow(1-tx ,2)/h ;
    }
    prod=prod*tp;
  }

  return prod;
}
struct QuarticKernel : public Worker
{
  // source matrix
  const RMatrix<double> input;

  double h;
  int p;
  // destination matrix
  RMatrix<double> output;

  // initialize with source and destination
  QuarticKernel(NumericMatrix input, NumericMatrix output, double h )
    : input(input),h(h), output(output) {}


  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++){
      for (std::size_t j = 0; j <= i; j++){
        RMatrix<double>::Row row1 = input.row(i);
        RMatrix<double>::Row row2 = input.row(j);
        // function to apply
        output(i,j) = Quartic_kernel(row1.begin(),row2.begin(),row1.end(),h);
        output(j,i) = output(i,j);
      }

    }
  }
};

NumericMatrix KernelQuartic_RcppParallel(NumericMatrix x, double h) {
  NumericMatrix rmat(x.nrow(), x.nrow());
  // declare the  instance
  QuarticKernel Kh_epan(x,rmat,h);
  // call parallelFor to start the work
  parallelFor(0, x.nrow(), Kh_epan);
  // return the computed matrix
  return rmat;
}


// triweight without boost
template <typename InputIterator1, typename InputIterator2>
inline double Triweight_kernel(InputIterator1 begin1,InputIterator2 begin2,
                               InputIterator1 end1, double h) {

  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  double prod=1.0;
  while(it1 != end1){
    double d1 = *it1++;
    double d2 = *it2++;
    double tx =pow((d1-d2)/h,2);
    double tp =0;
    if ( tx >= 1 ) {
      tp = 0;
    }else{
      tp =1.09375*pow(1-tx ,3)/h ;
    }
    prod=prod*tp;
  }

  return prod;
}
struct TriweightKernel : public Worker
{
  // source matrix
  const RMatrix<double> input;

  double h;
  int p;
  // destination matrix
  RMatrix<double> output;

  // initialize with source and destination
  TriweightKernel(NumericMatrix input, NumericMatrix output, double h )
    : input(input),h(h), output(output) {}


  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++){
      for (std::size_t j = 0; j <= i; j++){
        RMatrix<double>::Row row1 = input.row(i);
        RMatrix<double>::Row row2 = input.row(j);
        // function to apply
        output(i,j) = Triweight_kernel(row1.begin(),row2.begin(),row1.end(),h);
        output(j,i) = output(i,j);
      }

    }
  }
};

NumericMatrix KernelTriweight_RcppParallel(NumericMatrix x, double h) {
  NumericMatrix rmat(x.nrow(), x.nrow());
  // declare the  instance
  TriweightKernel Kh_epan(x,rmat,h);
  // call parallelFor to start the work
  parallelFor(0, x.nrow(), Kh_epan);
  // return the computed matrix
  return rmat;
}

// triangle without boost
template <typename InputIterator1, typename InputIterator2>
inline double Triangle_kernel(InputIterator1 begin1,InputIterator2 begin2,
                              InputIterator1 end1, double h) {

  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  double prod=1.0;
  while(it1 != end1){
    double d1 = *it1++;
    double d2 = *it2++;
    double tx =abs((d1-d2)/h);
    double tp =0;
    if ( tx >= 1 ) {
      tp = 0;
    }else{
      tp =(1-tx)/h ;
    }
    prod=prod*tp;
  }

  return prod;
}
struct TriangleKernel : public Worker
{
  // source matrix
  const RMatrix<double> input;

  double h;
  int p;
  // destination matrix
  RMatrix<double> output;

  // initialize with source and destination
  TriangleKernel(NumericMatrix input, NumericMatrix output, double h )
    : input(input),h(h), output(output) {}


  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++){
      for (std::size_t j = 0; j <= i; j++){
        RMatrix<double>::Row row1 = input.row(i);
        RMatrix<double>::Row row2 = input.row(j);
        // function to apply
        output(i,j) = Triangle_kernel(row1.begin(),row2.begin(),row1.end(),h);
        output(j,i) = output(i,j);
      }

    }
  }
};
NumericMatrix KernelTriangle_RcppParallel(NumericMatrix x, double h) {
  NumericMatrix rmat(x.nrow(), x.nrow());
  // declare the  instance
  TriangleKernel Kh_epan(x,rmat,h);
  // call parallelFor to start the work
  parallelFor(0, x.nrow(), Kh_epan);
  // return the computed matrix
  return rmat;
}


// cosine without boost
template <typename InputIterator1, typename InputIterator2>
inline double Cosine_kernel(InputIterator1 begin1,InputIterator2 begin2,
                            InputIterator1 end1, double h) {

  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  double prod=1.0;
  while(it1 != end1){
    double d1 = *it1++;
    double d2 = *it2++;
    double pi = pi = 3.14159265358979323846;
    double tx =abs((d1-d2)/h);
    double tp =0;
    if ( tx >= 1 ) {
      tp = 0;
    }else{
      tp =pi/4.0 * cos(pi*tx/2)/h ;
    }
    prod=prod*tp;
  }

  return prod;
}
struct CosineKernel : public Worker
{
  // source matrix
  const RMatrix<double> input;

  double h;
  int p;
  // destination matrix
  RMatrix<double> output;

  // initialize with source and destination
  CosineKernel(NumericMatrix input, NumericMatrix output, double h )
    : input(input),h(h), output(output) {}


  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++){
      for (std::size_t j = 0; j <= i; j++){
        RMatrix<double>::Row row1 = input.row(i);
        RMatrix<double>::Row row2 = input.row(j);
        // function to apply
        output(i,j) = Cosine_kernel(row1.begin(),row2.begin(),row1.end(),h);
        output(j,i) = output(i,j);
      }

    }
  }
};

NumericMatrix KernelCosine_RcppParallel(NumericMatrix x, double h) {
  NumericMatrix rmat(x.nrow(), x.nrow());
  // declare the  instance
  CosineKernel Kh_epan(x,rmat,h);
  // call parallelFor to start the work
  parallelFor(0, x.nrow(), Kh_epan);
  // return the computed matrix
  return rmat;
}

// [[Rcpp::export]]
NumericMatrix multi_kernelmatrix(NumericMatrix X, double h, std::string kernel){

  if(kernel=="gaussian"||kernel=="g"){
    return KernelGaussian_RcppParallel(X,h) ;
  }else if(kernel=="epanechnikov"||kernel=="e"){
    return KernelEpanechnikov_RcppParallel(X,h);
  }else if(kernel=="quartic"||kernel=="q"){
    return KernelQuartic_RcppParallel(X,h);
  } else if(kernel=="triweight"||kernel=="triw"){
    return KernelTriweight_RcppParallel(X,h);
  }  else if(kernel=="triangle"||kernel=="tria"){
    return KernelTriangle_RcppParallel(X,h);
  }  else if(kernel=="cosine"||kernel=="c"){
    return KernelCosine_RcppParallel(X,h);
  } else {
    Rcout << "Not included! \n";
    return -1;
  }

}
