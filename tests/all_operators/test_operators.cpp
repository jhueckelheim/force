#include <iostream>
#include <assert.h>
#include <random>
#include "force.hpp"

int main(int argc, char** argv) {
  int fixed_op = -1;
  int n_ops = 10;
  if(argc >= 2) {
    n_ops = atoi(argv[1]);
  }
  if(argc >= 3) {
    fixed_op = atoi(argv[2]);
  }
  std::cout<<"n_ops: "<<n_ops<<"; fixed_op: "<<fixed_op<<std::endl;
  // Test Constructors
  {
    freal<float> a, b(1.0), c(1.0,2.0);
    assert(a.value() == 0.0 && a.error() == 0.0 && a.corrected_value() == 0.0);
    assert(b.value() == 1.0 && b.error() == 0.0 && b.corrected_value() == 1.0);
    assert(c.value() == 1.0 && c.error() == 2.0 && c.corrected_value() == -1.0);
  }

  // Test Boolean Operators
  {
    freal<float> a(2.0,1.0), b(1.0), c(2.0,0.0);
    assert(a == b);
    assert(a != c);
    assert(a < c);
    assert(!(c < a));
    assert(c > a);
    assert(!(a > c));
    assert(a <= c);
    assert(!(a >= c));
    assert(c >= a);
    assert(!(c <= a));
    assert(a <= b);
    assert(a >= b);
  }
  // Test + - * / += -+ *= /= sin cos
  // We start with three random numbers in [-1,1], then randomly assign them to
  // A,B,C. Then we randomly pick (op) from [+ - * / sin cos += -= *= /=]
  // and perform "C = A (op) B". We do this 100 times. After each operation, we
  // divide C by a random value in [0,|C|] to avoid divergence.
  {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<float> rdist(-1,1);
    std::uniform_int_distribution<> dice6(0, 5);
    std::uniform_int_distribution<> dice10(0, 9);
    float old_as, old_bs, old_cs, new_as, new_bs, new_cs;
    double old_ad, old_bd, old_cd, new_ad, new_bd, new_cd;
    freal<float> old_af, old_bf, old_cf, new_af, new_bf, new_cf;
    new_as = rdist(generator);
    new_bs = rdist(generator);
    new_cs = rdist(generator);
    new_ad = (double) new_as;
    new_bd = (double) new_bs;
    new_cd = (double) new_cs;
    new_af = freal<float>(new_as);
    new_bf = freal<float>(new_bs);
    new_cf = freal<float>(new_cs);
    std::cout<<"starting with "<<new_as<<" "<<new_bs<<" "<<new_cs<<std::endl;
    for(int i=0; i<n_ops; i++) {
      int argpermute = dice6(generator);
      old_af = new_af; old_as = new_as; old_ad = new_ad;
      old_bf = new_bf; old_bs = new_bs; old_bd = new_bd;
      old_cf = new_cf; old_cs = new_cs; old_cd = new_cd;
      switch(argpermute) {
        case 0:
          new_af = old_af; new_as = old_as; new_ad = old_ad;
          new_bf = old_bf; new_bs = old_bs; new_bd = old_bd;
          new_cf = old_cf; new_cs = old_cs; new_cd = old_cd;
          break;
        case 1:
          new_af = old_af; new_as = old_as; new_ad = old_ad;
          new_bf = old_cf; new_bs = old_cs; new_bd = old_cd;
          new_cf = old_bf; new_cs = old_bs; new_cd = old_bd;
          break;
        case 2:
          new_af = old_bf; new_as = old_bs; new_ad = old_bd;
          new_bf = old_af; new_bs = old_as; new_bd = old_ad;
          new_cf = old_cf; new_cs = old_cs; new_cd = old_cd;
          break;
        case 3:
          new_af = old_bf; new_as = old_bs; new_ad = old_bd;
          new_bf = old_cf; new_bs = old_cs; new_bd = old_cd;
          new_cf = old_af; new_cs = old_as; new_cd = old_ad;
          break;
        case 4:
          new_af = old_cf; new_as = old_cs; new_ad = old_cd;
          new_bf = old_af; new_bs = old_as; new_bd = old_ad;
          new_cf = old_bf; new_cs = old_bs; new_cd = old_bd;
          break;
        case 5:
          new_af = old_cf; new_as = old_cs; new_ad = old_cd;
          new_bf = old_bf; new_bs = old_bs; new_bd = old_bd;
          new_cf = old_af; new_cs = old_as; new_cd = old_ad;
          break;
      }
      int op = dice10(generator);
      if(fixed_op != -1) {
        op = fixed_op;
      }
      switch(op) {
        case 0:
          new_cf = new_af + new_bf; new_cs = new_as + new_bs; new_cd = new_ad + new_bd;
          break;
        case 1:
          new_cf = new_af - new_bf; new_cs = new_as - new_bs; new_cd = new_ad - new_bd;
          break;
        case 2:
          new_cf = new_af * new_bf; new_cs = new_as * new_bs; new_cd = new_ad * new_bd;
          break;
        case 3:
          if(new_bd != 0 && new_bs != 0 && new_bf.value() != 0) {
            new_cf = new_af / new_bf; new_cs = new_as / new_bs; new_cd = new_ad / new_bd;
          }
          break;
        case 4:
          new_cf += new_bf; new_cs += new_bs; new_cd += new_bd;
          break;
        case 5:
          new_cf -= new_bf; new_cs -= new_bs; new_cd -= new_bd;
          break;
        case 6:
          new_cf *= new_bf; new_cs *= new_bs; new_cd *= new_bd;
          break;
        case 7:
          if(new_bd != 0 && new_bs != 0 && new_bf.value() != 0) {
            new_cf /= new_bf; new_cs /= new_bs; new_cd /= new_bd;
          }
          break;
        case 8:
          new_cf = sin(new_bf); new_cs = sin(new_bs); new_cd = sin(new_bd);
          break;
        case 9:
          new_cf = cos(new_bf); new_cs = cos(new_bs); new_cd = cos(new_bd);
          break;
      }
    }
    std::cout<<"float : "<<new_as<<" "<<new_bs<<" "<<new_cs<<std::endl;
    std::cout<<"double: "<<new_ad<<" "<<new_bd<<" "<<new_cd<<std::endl;
    std::cout<<"force : "<<new_af<<" "<<new_bf<<" "<<new_cf<<std::endl;
  }
  return 0;
}
