#ifndef force_hpp
#define force_hpp

#include <iostream>
#include <math.h>
#ifdef MPFR
#include "mpfrcpp_tpl.h"
#endif

template<typename T>
class freal{
  private:
    T val, err;
    #ifdef _CLASSIC_CENA_ADDITION
    static T addition_error(T a, T b, T x) {
        if (fabs(a) >= fabs(b)) {
            return (x - a) - b;
        } else {
            return (x - b) - a;
        }
    }
    #else
    static T addition_error(T a, T b, T x) {
        // addition error without branching, as in Ogita et al: "ACCURATE SUM AND DOT PRODUCT"
        T corr = x - a;
        return ((x-corr)-a)+(corr-b);
    }
    #endif

    static T multiplication_error(T a, T b, T x, T m) {
      T au, al, bu, bl;
      au = (a - a*m) + a*m;
      al = a - au;
      bu = (b - b*m) + b*m;
      bl = b - bu;
      return x - au*bu - (au*bl + al*bu) - al*bl;
    }

    static float multiplication_error(float a, float b, float x) {
      return multiplication_error(a,b,x,(float)4097); //pow(2,mantissalength/2) + 1;
    }

    static double multiplication_error(double a, double b, double x) {
      return multiplication_error(a,b,x,(double)67108865); //pow(2,mantissalength/2) + 1;
    }

    static long double multiplication_error(long double a, long double b, long double x) {
      return multiplication_error(a,b,x,(long double)4294967297.0L); //pow(2,mantissalength/2) + 1;
    }

#ifdef MPFR
    template<unsigned int prec>
    static mpfrcpp<prec> multiplication_error(mpfrcpp<prec> a, mpfrcpp<prec> b, mpfrcpp<prec> x) {
      return multiplication_error(a,b,x,pow(mpfrcpp<prec>(2.0),prec/2)+mpfrcpp<prec>(1.0)); //pow(2,mantissalength/2) + 1;
    }
#endif

    static T division_error(T a, T b, T x) {
      return (x*b - a - multiplication_error(x, b, x*b)) / b;
    }

    static T sqrt_error(T a, T x) {
      return (x*x - a + multiplication_error(x, x, x*x)) / (((T)2.0)*x);
    }

  public:
    freal<T>(T value, T error) : val(value), err(error) { }
    freal<T>(T value) : val(value), err((T)0.0) { }
    freal<T>() : val((T)0.0), err((T)0.0) { }
    T value() const { return this->val; }
    T error() const { return this->err; }
    T corrected_value() const { return this->val-this->err; }
    operator T() const { return corrected_value(); }

    freal<T>& operator=(const freal<T> &g1){
      this->val = g1.val;
      this->err = g1.err;
      return *this;
    }
    freal<T>& operator=(const T &g1){
      this->val = g1;
      this->err = (T)0.0;
      return *this;
    }
    void operator+=(const freal<T> rhs) {
      T value = this->val + rhs.val;
      T localerror = addition_error(this->val,rhs.val,value);
      this->val = value;
      this->err += rhs.err + localerror;
    }
    void operator-=(const freal<T> rhs) {
      T value = this->val - rhs.val;
      T localerror = addition_error(this->val,-rhs.val,value);
      this->val = value;
      this->err -= rhs.err + localerror;
    }
    void operator*=(const freal<T> rhs) {
      T value = this->val * rhs.val;
      T localerr = multiplication_error(this->val,rhs.val,value);
      this->err = this->val*rhs.err+this->err*rhs.val + localerr;
      this->val = value;
    }
    void operator/=(const freal<T> rhs) {
      T recip = 1.0 / rhs.val;
      T newval = this->val / rhs.val;
      T localerr = division_error(this->val,rhs.val,newval);
      this->val = newval;
      this->err = recip*this->err - recip*newval*rhs.err + localerr;
    }
    static void ompReduce(freal<T> &omp_out, freal<T> &omp_in) {
      omp_out += omp_in;
    }

  template <typename T1>
  friend freal<T1> operator+(const freal<T1> &g1,const freal<T1> &g2);
  template <typename T1>
  friend freal<T1> operator-(const freal<T1> &g1,const freal<T1> &g2);
  template <typename T1>
  friend freal<T1> operator-(const freal<T1> &g1);
  template <typename T1>
  friend freal<T1> operator*(const freal<T1> &g1,const freal<T1> &g2);
  template <typename T1>
  friend freal<T1> operator/(const freal<T1> &g1,const freal<T1> &g2);
  template <typename T1>
  friend bool operator==(const freal<T1> &g1,const freal<T1> &g2);
  template <typename T1>
  friend bool operator!=(const freal<T1> &g1,const freal<T1> &g2);
  template <typename T1>
  friend bool operator>(const freal<T1> &g1,const freal<T1> &g2);
  template <typename T1>
  friend bool operator<(const freal<T1> &g1,const freal<T1> &g2);
  template <typename T1>
  friend bool operator>=(const freal<T1> &g1,const freal<T1> &g2);
  template <typename T1>
  friend bool operator<=(const freal<T1> &g1,const freal<T1> &g2);
  template <typename T1>
  friend freal<T1> sin(const freal<T1> &);
  template <typename T1>
  friend freal<T1> cos(const freal<T1> &);
  template <typename T1>
  friend freal<T1> sqrt(const freal<T1> &);
  template <typename T1>
  friend std::ostream& operator<<(std::ostream&, const freal<T1>&);
};

template<typename T>
freal<T> operator+(const freal<T> &g1,const freal<T> &g2){
  T value = g1.val + g2.val;
  T localerror = freal<T>::addition_error(g1.val,g2.val,value);
  T error = g1.err + g2.err + localerror;
  return freal<T>(value,error);
}

template<typename T>
freal<T> operator-(const freal<T> &g1,const freal<T> &g2){
  T value = g1.val - g2.val;
  T localerror = freal<T>::addition_error(g1.val,-g2.val,value);
  T error = g1.err - g2.err + localerror;
  return freal<T>(value,error);
}

template<typename T>
freal<T> operator-(const freal<T> &g1){
   return freal<T>(-g1.val,g1.err);
}

template<typename T>
freal<T> operator*(const freal<T> &g1,const freal<T> &g2){
   T value = g1.val * g2.val;
   T localerr = freal<T>::multiplication_error(g1.val,g2.val,value);
   T error = g1.val*g2.err+g1.err*g2.val + localerr;
   return freal<T>(value,error);
}

template<typename T>
freal<T> operator/(const freal<T> &g1,const freal<T> &g2){
   T recip,newval;
   recip = (T)1.0 / g2.val;
   newval = g1.val / g2.val;
   T localerr = freal<T>::division_error(g1.val,g2.val,newval);
   return freal<T>(newval,recip*g1.err - recip*newval*g2.err + localerr);
}

template<typename T>
bool operator==(const freal<T> &g1,const freal<T> &g2){
   return g1.corrected_value() == g2.corrected_value();
}

template<typename T>
bool operator!=(const freal<T> &g1,const freal<T> &g2){
   return g1.corrected_value() != g2.corrected_value();
}

template<typename T>
bool operator>(const freal<T> &g1,const freal<T> &g2){
   return g1.corrected_value() > g2.corrected_value();
}

template<typename T>
bool operator<(const freal<T> &g1,const freal<T> &g2){
   return g1.corrected_value() < g2.corrected_value();
}

template<typename T>
bool operator>=(const freal<T> &g1,const freal<T> &g2){
   return g1.corrected_value() >= g2.corrected_value();
}

template<typename T>
bool operator<=(const freal<T> &g1,const freal<T> &g2){
   return g1.corrected_value() <= g2.corrected_value();
}

template<typename T>
freal<T> sin(const freal<T> &g1){
  return freal<T>(sin(g1.val),cos(g1.val)*g1.err);
}

template<typename T>
freal<T> cos(const freal<T> &g1){
  return freal<T>(cos(g1.val),-sin(g1.val)*g1.err);
}

template<typename T>
freal<T> fabs(const freal<T> &g1){
  return g1>=freal<float>(0.0)?g1:-g1;
}

template<typename T>
freal<T> sqrt(const freal<T> &g1){
   T newval,two;
   two = 2.0;
   newval = sqrt(g1.val);
   T localerr = freal<T>::sqrt_error(g1.val,newval);
   return freal<T>(newval,g1.err/(two*newval) + localerr);
}

template<typename T>
std::ostream& operator<<(std::ostream &ost, const freal<T> &ad){
  ost << "[" << ad.value() << ", " << ad.error() << ", " << (ad.corrected_value()) << "]";
  return ost;
}

#pragma omp declare reduction(+ : freal<float> : freal<float>::ompReduce(omp_out, omp_in))
#pragma omp declare reduction(+ : freal<double> : freal<double>::ompReduce(omp_out, omp_in)) 
#pragma omp declare reduction(+ : freal<long double> : freal<long double>::ompReduce(omp_out, omp_in)) 

#endif
