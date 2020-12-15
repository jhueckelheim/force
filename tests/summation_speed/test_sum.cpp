#include <iostream>
#include <iomanip>
#include "force.hpp"
#include "quadhelper.hpp"
#include <random>
#include <chrono>

template<typename T>
class KahanReal {
  public:
    T sum, cor;
    KahanReal<T>() : sum(0), cor(0) {}
    static void ompAdd(KahanReal<T> &omp_out, KahanReal<T> &omp_in) {
      T y = omp_in.sum - (omp_out.cor + omp_in.cor);
      T t = omp_out.sum + y;
      omp_out.cor = (t - omp_out.sum) - y;
      omp_out.sum = t;
    }
};
#pragma omp declare reduction(+ : KahanReal<float> : KahanReal<float>::ompAdd(omp_out, omp_in))
#pragma omp declare reduction(+ : KahanReal<double> : KahanReal<double>::ompAdd(omp_out, omp_in))
#pragma omp declare reduction(+ : KahanReal<long double> : KahanReal<long double>::ompAdd(omp_out, omp_in))

template<typename T>
T benchmark(double* dval, int n, double& time) {
  auto val = new T[n];
  for(int i=0; i<n; i++) {
    val[i] = (T)dval[i];
  }
  auto t_start = std::chrono::high_resolution_clock::now();
  T sum = 0;
  #pragma omp parallel for reduction(+:sum)
  for(int i=0; i<n; i++) {
    sum += val[i];
  }
  auto t_end = std::chrono::high_resolution_clock::now();
  time = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  delete[] val;
  return sum;
}

template<typename T>
T benchmark_kahan(double* dval, int n, double& time) {
  auto val = new T[n];
  for(int i=0; i<n; i++) {
    val[i] = (T)dval[i];
  }
  auto t_start = std::chrono::high_resolution_clock::now();
  KahanReal<T> sum;
  #pragma omp parallel for reduction(+:sum)
  for(int i=0; i<n; i++) {
    T y = val[i] - sum.cor;
    T t = sum.sum + y;
    sum.cor = (t - sum.sum) - y;
    sum.sum = t;
  }
  auto t_end = std::chrono::high_resolution_clock::now();
  time = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  delete[] val;
  return sum.sum;
}

int main(int argc, char** argv) {
  // Read command line arguments and set up random input data
  if(argc < 3) {
    std::cout<<"Usage: exec_test_sum <numvals> <distribution>"<<std::endl;
    exit(-1);
  }
  int n = atoi(argv[1]);
  float dist = atof(argv[2]);
  auto dval = new double[n];
  std::default_random_engine generator;
  std::exponential_distribution<double> distribution(dist);
  double elapsed_time;
  for(int i=0; i<n; i++) {
    dval[i] = distribution(generator);
  }
  // Increase output precision for floating point numbers
  std::cout<<std::setprecision(36);
  // Warmup run to make timings consistent
  {
    float sum = benchmark<float>(dval, n, elapsed_time);
    std::cout<<"warmup "<<sum<<" time "<<elapsed_time<<std::endl;
  }

  // Benchmark single precision
  {
    float sum = benchmark<float>(dval, n, elapsed_time);
    std::cout<<"single_sum "<<sum<<" time "<<elapsed_time<<std::endl;
  }
  // Benchmark double precision
  {
    double sum = benchmark<double>(dval, n, elapsed_time);
    std::cout<<"double_sum "<<sum<<" time "<<elapsed_time<<std::endl;
  }
  // Benchmark long double precision
  {
    long double sum = benchmark<long double>(dval, n, elapsed_time);
    std::cout<<"longdouble_sum "<<sum<<" time "<<elapsed_time<<std::endl;
  }
  // Benchmark quad precision
  {
    __float128 sum = benchmark<__float128>(dval, n, elapsed_time);
    std::cout<<"quad_sum "<<sum<<" time "<<elapsed_time<<std::endl;
  }
  // Benchmark single precision with CENA
  {
    freal<float> sum = benchmark<freal<float>>(dval, n, elapsed_time);
    std::cout<<"csingle_sum "<<sum<<" time "<<elapsed_time<<std::endl;
  }
  // Benchmark double precision with CENA
  {
    freal<double> sum = benchmark<freal<double>>(dval, n, elapsed_time);
    std::cout<<"cdouble_sum "<<sum<<" time "<<elapsed_time<<std::endl;
  }
  // Benchmark long double precision with CENA
  {
    freal<long double> sum = benchmark<freal<long double>>(dval, n, elapsed_time);
    std::cout<<"clongdouble_sum "<<sum<<" time "<<elapsed_time<<std::endl;
  }
  // Benchmark single precision with Kahan
  {
    float sum = benchmark_kahan<float>(dval, n, elapsed_time);
    std::cout<<"ksingle_sum "<<sum<<" time "<<elapsed_time<<std::endl;
  }
  // Benchmark double precision with Kahan
  {
    double sum = benchmark_kahan<double>(dval, n, elapsed_time);
    std::cout<<"kdouble_sum "<<sum<<" time "<<elapsed_time<<std::endl;
  }
  // Benchmark long double precision with Kahan
  {
    long double sum = benchmark_kahan<long double>(dval, n, elapsed_time);
    std::cout<<"klongdouble_sum "<<sum<<" time "<<elapsed_time<<std::endl;
  }
  //// Benchmark single precision Kahan summation
  //{
  //  auto t_start = std::chrono::high_resolution_clock::now();
  //  float sum = 0;
  //  float cor = 0;
  //  for(int i=0; i<n; i++) {
  //    float y = sval[i] - cor;
  //    float t = sum + y;
  //    cor = (t - sum) - y;
  //    sum = t;
  //  }
  //  auto t_end = std::chrono::high_resolution_clock::now();
  //  double elapsed = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  //  std::cout<<"kahan_serl "<<sum<<" time "<<elapsed<<std::endl;
  //}
  //// Benchmark single precision parallel Kahan summation
  //{
  //  auto t_start = std::chrono::high_resolution_clock::now();
  //  KahanReal<float> sum;
  //  #pragma omp parallel for reduction(+:sum)
  //  for(int i=0; i<n; i++) {
  //    float y = sval[i] - sum.cor;
  //    float t = sum.sum + y;
  //    sum.cor = (t - sum.sum) - y;
  //    sum.sum = t;
  //  }
  //  auto t_end = std::chrono::high_resolution_clock::now();
  //  double elapsed = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  //  std::cout<<"kahan_prll "<<sum.sum<<" time "<<elapsed<<std::endl;
  //}

  // Clean up and exit.
  delete[] dval;
  return 0;
}

