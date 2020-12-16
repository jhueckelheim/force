#ifndef mpfrcpp_h
#define mpfrcpp_h

#include <iostream>
#include <gmp.h>
#include <mpfr.h>
#include <quadmath.h>
#include <stdlib.h>
#include <stdio.h>
template <unsigned int MPFRPREC>
class mpfrcpp{
public:
   mpfr_t value;
   mpfrcpp() {
      mpfr_init2(value,MPFRPREC);
   }
   mpfrcpp(const float v) {
      mpfr_init2(value,MPFRPREC);
      mpfr_set_flt(value,v,MPFR_RNDN);
   }
   mpfrcpp(const double v) {
      mpfr_init2(value,MPFRPREC);
      mpfr_set_d(value,v,MPFR_RNDN);
   }
   mpfrcpp(const long double v) {
      mpfr_init2(value,MPFRPREC);
      mpfr_set_ld(value,v,MPFR_RNDN);
   }
   mpfrcpp(const __float128 v) {
      char buf[128];
      quadmath_snprintf(buf, sizeof(buf), "%.36Qg", v);
      mpfr_init2(value,MPFRPREC);
      mpfr_set_str(value,buf,10,MPFR_RNDN);
   }
   mpfrcpp(const mpfrcpp &v) {
      mpfr_init2(value,MPFRPREC);
      mpfr_set(value,v.value,MPFR_RNDN);
   }
   mpfrcpp(const mpfr_t &v) {
      mpfr_init2(value,MPFRPREC);
      mpfr_set(value,v,MPFR_RNDN);
   }
   mpfrcpp(const char* v) {
      mpfr_init2(value,MPFRPREC);
      mpfr_set_str(value,v,10,MPFR_RNDN);
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

template<unsigned int fromprec, unsigned int toprec>
mpfrcpp<toprec> convert(const mpfrcpp<fromprec> &from) {
   mpfr_t to;
   mpfr_init2(to,toprec);
   mpfr_set(to,from.value,MPFR_RNDN);
   return to;
}

template<unsigned int MPFRPREC>
mpfrcpp<MPFRPREC> operator/(const mpfrcpp<MPFRPREC> &g1,const mpfrcpp<MPFRPREC> &g2){
   mpfrcpp<MPFRPREC> res;
   mpfr_div(res.value,g1.value,g2.value,MPFR_RNDN);
   return res;
}
template<unsigned int MPFRPREC>
std::ostream& operator<<(std::ostream &ost, const mpfrcpp<MPFRPREC> &ad){
   char* abc = NULL;
   mpfr_exp_t i;
   if(ad >= mpfrcpp<200>(0.0)) {
     abc = mpfr_get_str (NULL, &i, 10, 0, ad.value, MPFR_RNDN);
     ost << "0." << abc << "e" << i;
   }
   else {
     abc = mpfr_get_str (NULL, &i, 10, 0, (-ad).value, MPFR_RNDN);
     ost << "-0." << abc << "e" << i;
   }
   mpfr_free_str(abc);
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

