#include <iostream>
#include <iomanip>
#include "force.hpp"
#include "mpfrcpp_tpl.h"
#include "quadhelper.hpp"
#include "randomhelper.hpp"
#include <random>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <unistd.h>
#include <random>
#include <iostream>

template<typename T>
T addup(std::vector<__float128> qval, int n) {
  auto val = new T[n];
  for(int i=0; i<n; i++) {
    val[i] = (T)qval[i];
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
T addup_kahan(std::vector<__float128> qval, int n) {
  auto val = new T[n];
  for(int i=0; i<n; i++) {
    val[i] = (T)qval[i];
  }
  T sum, cor;
  sum = 0; cor = 0;
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
  std::vector<__float128> qval;
  // Here we generate "random" numbers from an exponential distribution.
  // Note that we do not seed the random engine, and hence the sequence
  // of generated numbers will always be the same between runs of this
  // program. We use a custom-made generator for high-resolution random
  // numbers where all mantissa bits are uniformly random, and the
  // exponent and sign bits are uniformly random from an allowable range.
  double elapsed_time;
  for(int i=0; i<n; i++) {
    qval.push_back(rand_quad());
  }
  // Now we shuffle the vector with numbers. This time we use a seeded
  // random number generator, so that we will add the numbers in a
  // different order each time this program is run. This seed is a bit
  // funky, since we run many benchmarks back-to-back and a seed based
  // only on a low-resolution clock would sometimes remain unchanged
  // between two runs.
  // Sadly, using a C++11 random_device still seems to rely on the time
  // to seed on some machines, so we fall back to good old srand here.
  srand((time(NULL) & 0xFFFF) | (getpid() << 16));
  std::random_shuffle ( qval.begin(), qval.end() );
  // Increase output precision for floating point numbers
  std::cout<<std::setprecision(36);

  // Compute high-precision reference result
  mpfrcpp<200> rsum = addup<mpfrcpp<200>>(qval, n);
  std::cout<<"reference_sum "<<rsum<<std::endl;
  // Benchmark single precision
  {
    float sum = addup<float>(qval, n);
    std::cout<<"single_sum "<<sum<<" err "<<rsum-mpfrcpp<200>(sum)<<std::endl;
  }
  // Benchmark double precision
  {
    double sum = addup<double>(qval, n);
    std::cout<<"double_sum "<<sum<<" err "<<rsum-mpfrcpp<200>(sum)<<std::endl;
  }
  // Benchmark long double precision
  {
    long double sum = addup<long double>(qval, n);
    std::cout<<"longdouble_sum "<<sum<<" err "<<rsum-mpfrcpp<200>(sum)<<std::endl;
  }
  // Benchmark quad precision
  {
    __float128 sum = addup<__float128>(qval, n);
    std::cout<<"quad_sum "<<sum<<" err "<<rsum-mpfrcpp<200>(sum)<<std::endl;
  }
  // Benchmark single precision with CENA
  {
    freal<float> sum = addup<freal<float>>(qval, n);
    std::cout<<"csingle_sum "<<(float)sum<<" err "<<rsum-mpfrcpp<200>(sum)<<std::endl;
  }
  // Benchmark double precision with CENA
  {
    freal<double> sum = addup<freal<double>>(qval, n);
    std::cout<<"cdouble_sum "<<(double)sum<<" err "<<rsum-mpfrcpp<200>(sum)<<std::endl;
  }
  // Benchmark long double precision with CENA
  {
    freal<long double> sum = addup<freal<long double>>(qval, n);
    std::cout<<"clongdouble_sum "<<(long double)sum<<" err "<<rsum-mpfrcpp<200>(sum)<<std::endl;
  }
  // Benchmark single precision with Kahan
  {
    float sum = addup_kahan<float>(qval, n);
    std::cout<<"ksingle_sum "<<sum<<" err "<<rsum-mpfrcpp<200>(sum)<<std::endl;
  }
  // Benchmark double precision with Kahan
  {
    double sum = addup_kahan<double>(qval, n);
    std::cout<<"kdouble_sum "<<sum<<" err "<<rsum-mpfrcpp<200>(sum)<<std::endl;
  }
  // Benchmark long double precision with Kahan
  {
    long double sum = addup_kahan<long double>(qval, n);
    std::cout<<"klongdouble_sum "<<sum<<" err "<<rsum-mpfrcpp<200>(sum)<<std::endl;
  }

  return 0;
}

