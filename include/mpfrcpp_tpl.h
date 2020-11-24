#ifndef mpfrcpp_h
#define mpfrcpp_h

#include <iostream>
#include <gmp.h>
#include <mpfr.h>
template <unsigned int MPFRPREC>
class mpfrcpp{
public:
   mpfr_t value;
   mpfrcpp() {
      mpfr_init2(value,MPFRPREC);
   }
   mpfrcpp(const double v) {
      mpfr_init2(value,MPFRPREC);
      mpfr_set_d(value,v,MPFR_RNDN);
   }
   mpfrcpp(const mpfrcpp &v) {
      mpfr_init2(value,MPFRPREC);
      mpfr_set(value,v.value,MPFR_RNDN);
   }
   mpfrcpp(const mpfr_t &v) {
      mpfr_init2(value,MPFRPREC);
      mpfr_set(value,v,MPFR_RNDN);
   }
   ~mpfrcpp() {
      mpfr_clear(value);
   }

   mpfrcpp& operator=(const mpfrcpp &g1){
      mpfr_set(value,g1.value,MPFR_RNDN);
      return *this;
   }
   mpfrcpp& operator=(const double &g1){
      mpfr_set_d(value,g1,MPFR_RNDN);
      return *this;
   }
};

template<unsigned int MPFRPREC>
mpfrcpp<MPFRPREC> operator/(const mpfrcpp<MPFRPREC> &g1,const mpfrcpp<MPFRPREC> &g2){
   mpfrcpp<MPFRPREC> res;
   mpfr_div(res.value,g1.value,g2.value,MPFR_RNDN);
   return res;
}
template<unsigned int MPFRPREC>
std::ostream& operator<<(std::ostream &ost, const mpfrcpp<MPFRPREC> &ad){
   ost << mpfr_get_d(ad.value,MPFR_RNDN);
   return ost;
}
template<unsigned int MPFRPREC>
mpfrcpp<MPFRPREC> operator+(const mpfrcpp<MPFRPREC> &g1,const mpfrcpp<MPFRPREC> &g2){
   mpfrcpp<MPFRPREC> res;
   mpfr_add(res.value,g1.value,g2.value,MPFR_RNDN);
   return res;
}
template<unsigned int MPFRPREC>
mpfrcpp<MPFRPREC> operator-(const mpfrcpp<MPFRPREC> &g1,const mpfrcpp<MPFRPREC> &g2){
   mpfrcpp<MPFRPREC> res;
   mpfr_sub(res.value,g1.value,g2.value,MPFR_RNDN);
   return res;
}
template<unsigned int MPFRPREC>
mpfrcpp<MPFRPREC> operator-(const mpfrcpp<MPFRPREC> &g1){
   mpfrcpp<MPFRPREC> res;
   mpfr_neg(res.value,g1.value,MPFR_RNDN);
   return res;
}
template<unsigned int MPFRPREC>
mpfrcpp<MPFRPREC> operator*(const mpfrcpp<MPFRPREC> &g1,const mpfrcpp<MPFRPREC> &g2){
   mpfrcpp<MPFRPREC> res;
   mpfr_mul(res.value,g1.value,g2.value,MPFR_RNDN);
   return res;
}
template<unsigned int MPFRPREC>
mpfrcpp<MPFRPREC> fabs(const mpfrcpp<MPFRPREC> &g1){
   mpfrcpp<MPFRPREC> res;
   mpfr_abs(res.value,g1.value,MPFR_RNDN);
   return res;
}
template<unsigned int MPFRPREC>
mpfrcpp<MPFRPREC> pow(const mpfrcpp<MPFRPREC> &g1, double expd){
   mpfrcpp<MPFRPREC> res;
   mpfrcpp<MPFRPREC> exp = expd;
   mpfr_pow(res.value,g1.value,exp.value,MPFR_RNDN);
   return res;
}
template<unsigned int MPFRPREC>
mpfrcpp<MPFRPREC> sqrt(const mpfrcpp<MPFRPREC> &g1){
   mpfrcpp<MPFRPREC> res;
   mpfr_sqrt(res.value,g1.value,MPFR_RNDN);
   return res;
}
template<unsigned int MPFRPREC>
int operator>=(const mpfrcpp<MPFRPREC> &g1,const mpfrcpp<MPFRPREC> &g2){
   return mpfr_greaterequal_p(g1.value,g2.value);
}
template<unsigned int MPFRPREC>
int operator>(const mpfrcpp<MPFRPREC> &g1,const mpfrcpp<MPFRPREC> &g2){
   return mpfr_greater_p(g1.value,g2.value);
}
template<unsigned int MPFRPREC>
int operator<(const mpfrcpp<MPFRPREC> &g1,const mpfrcpp<MPFRPREC> &g2){
   return mpfr_greater_p(g2.value,g1.value);
}

#endif

