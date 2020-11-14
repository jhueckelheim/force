#include <iostream>
#include "force.hpp"
#include <random>

std::default_random_engine generator;

double exprand(float range) {
  std::exponential_distribution<double> distribution(range);
  return distribution(generator);
}

int main(int argc, char** argv) {
  if(argc < 3) {
    std::cout<<"Usage: exec_test_sum <numvals> <valrange>"<<std::endl;
    exit(-1);
  }
  int n = atoi(argv[1]);
  float range = atof(argv[2]);
  float ssum = 0;
  double dsum = 0;
  freal<float> fsum = 0;
  for(int i=0; i<n; i++) {
    double dval = exprand(range);
    float  sval = (float)dval;
    freal<float> fval = freal<float>(dval);
    ssum = ssum + sval;
    dsum = dsum + dval;
    fsum = fsum + fval;
  }
  std::cout<<"ssum: "<<ssum<<std::endl;
  std::cout<<"dsum: "<<dsum<<std::endl;
  std::cout<<"fsum: "<<fsum<<std::endl;
  return 0;
}


