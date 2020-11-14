#include <iostream>
#include "force.hpp"
#include <random>
#include <chrono>

int main(int argc, char** argv) {
  // Read command line arguments and set up random input data
  if(argc < 3) {
    std::cout<<"Usage: exec_test_sum <numvals> <distribution>"<<std::endl;
    exit(-1);
  }
  int n = atoi(argv[1]);
  float dist = atof(argv[2]);
  auto dval = new double[n];
  auto sval = new float[n];
  auto fval = new freal<float>[n];
  std::default_random_engine generator;
  std::exponential_distribution<double> distribution(dist);
  for(int i=0; i<n; i++) {
    dval[i] = distribution(generator);
    sval[i] = (float)dval[i];
    fval[i] = freal<float>(dval[i]);
  }

  // Benchmark single precision
  {
    auto t_start = std::chrono::high_resolution_clock::now();
    float sum = 0;
    for(int i=0; i<n; i++) {
      sum += sval[i];
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double, std::milli>(t_end-t_start).count();
    std::cout<<"single_sum "<<sum<<" time "<<elapsed<<std::endl;
  }
  // Benchmark double precision
  {
    auto t_start = std::chrono::high_resolution_clock::now();
    double sum = 0;
    for(int i=0; i<n; i++) {
      sum += dval[i];
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double, std::milli>(t_end-t_start).count();
    std::cout<<"double_sum "<<sum<<" time "<<elapsed<<std::endl;
  }
  // Benchmark single precision with Forward CENA
  {
    auto t_start = std::chrono::high_resolution_clock::now();
    freal<float> sum = 0;
    for(int i=0; i<n; i++) {
      sum += fval[i];
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double, std::milli>(t_end-t_start).count();
    std::cout<<"force_sum  "<<sum<<" time "<<elapsed<<std::endl;
  }
  // Benchmark single precision Kahan summation
  {
    auto t_start = std::chrono::high_resolution_clock::now();
    float sum = 0;
    float cor = 0;
    for(int i=0; i<n; i++) {
      float y = sval[i] - cor;
      float t = sum + y;
      cor = (t - sum) - y;
      sum = t;
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double, std::milli>(t_end-t_start).count();
    std::cout<<"kahan_sum "<<sum<<" time "<<elapsed<<std::endl;
  }

  // Clean up and exit.
  delete[] fval;
  delete[] sval;
  delete[] dval;
  return 0;
}

