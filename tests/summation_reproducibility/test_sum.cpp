#include <iostream>
#include <iomanip>
#include "force.hpp"
#include "mpfrcpp_tpl.h"
#include "quadhelper.hpp"
#include <random>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <unistd.h>

template<typename T>
T addup(std::vector<long double> dval, int n) {
  auto val = new T[n];
  for(int i=0; i<n; i++) {
    val[i] = (T)dval[i];
  }
  T sum;
  sum = 0;
  for(int i=0; i<n; i++) {
    sum = sum + val[i];
  }
  delete[] val;
  return sum;
}

template<typename T>
T addup_kahan(std::vector<long double> dval, int n) {
  auto val = new T[n];
  for(int i=0; i<n; i++) {
    val[i] = (T)dval[i];
  }
  T sum, cor;
  for(int i=0; i<n; i++) {
    T y = val[i] - cor;
    T t = sum + y;
    cor = (t - sum) - y;
    sum = t;
  }
  delete[] val;
  return sum;
}

int main(int argc, char** argv) {
  // Read command line arguments and set up random input data
  if(argc < 3) {
    std::cout<<"Usage: exec_test_sum <numvals> <distribution>"<<std::endl;
    exit(-1);
  }
  int n = atoi(argv[1]);
  float dist = atof(argv[2]);
  std::vector<long double> dval;
  // Here we generate "random" numbers from an exponential distribution.
  // Note that we do not seed the random engine, and hence the sequence
  // of generated numbers will always be the same between runs of this
  // program. We use an exponential distribution to get numbers from a
  // wide range of sizes.
  std::default_random_engine generator;
  std::exponential_distribution<long double> distribution(dist);
  double elapsed_time;
  for(int i=0; i<n; i++) {
    dval.push_back(distribution(generator));
  }
  // Now we shuffle the vector with numbers. This time we use a seeded
  // random number generator, so that we will add the numbers in a
  // different order each time this program is run. This seed is a bit
  // funky, since we run many benchmarks back-to-back and a seed based
  // only on a low-resolution clock would sometimes remain unchanged
  // between two runs.
  srand((time(NULL) & 0xFFFF) | (getpid() << 16));
  std::random_shuffle ( dval.begin(), dval.end() );
  // Increase output precision for floating point numbers
  std::cout<<std::setprecision(36);

  // Benchmark single precision
  {
    float sum = addup<float>(dval, n);
    std::cout<<"single_sum "<<sum<<std::endl;
  }
  // Benchmark double precision
  {
    double sum = addup<double>(dval, n);
    std::cout<<"double_sum "<<sum<<std::endl;
  }
  // Benchmark long double precision
  {
    long double sum = addup<long double>(dval, n);
    std::cout<<"longdouble_sum "<<sum<<std::endl;
  }
  // Benchmark quad precision
  {
    __float128 sum = addup<__float128>(dval, n);
    std::cout<<"quad_sum "<<sum<<std::endl;
  }
  // Benchmark single precision with CENA
  {
    freal<float> sum = addup<freal<float>>(dval, n);
    std::cout<<"csingle_sum "<<sum<<std::endl;
  }
  // Benchmark double precision with CENA
  {
    freal<double> sum = addup<freal<double>>(dval, n);
    std::cout<<"cdouble_sum "<<sum<<std::endl;
  }
  // Benchmark long double precision with CENA
  {
    freal<long double> sum = addup<freal<long double>>(dval, n);
    std::cout<<"clongdouble_sum "<<sum<<std::endl;
  }
  // Benchmark single precision with Kahan
  {
    float sum = addup_kahan<float>(dval, n);
    std::cout<<"ksingle_sum "<<sum<<std::endl;
  }
  // Benchmark double precision with Kahan
  {
    double sum = addup_kahan<double>(dval, n);
    std::cout<<"kdouble_sum "<<sum<<std::endl;
  }
  // Benchmark long double precision with Kahan
  {
    long double sum = addup_kahan<long double>(dval, n);
    std::cout<<"klongdouble_sum "<<sum<<std::endl;
  }
  // Benchmark long double precision with Kahan
  {
    mpfrcpp<200> sum = addup<mpfrcpp<200>>(dval, n);
    std::cout<<"reference_sum "<<sum<<std::endl;
  }

  return 0;
}

