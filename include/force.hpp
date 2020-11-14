#include <iostream>
#include <math.h>

template<typename T>
class freal{
  private:
    T value, error;
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

    static T multiplication_error(T a, T b, T x, int m) {
      T au, al, bu, bl;
      au = (a - a*m) + a*m;
      al = a - au;
      bu = (b - b*m) + b*m;
      bl = b - bu;
      return x - au*bu - (au*bl + al*bu) - al*bl;
    }

    static float multiplication_error(float a, float b, float x) {
      return multiplication_error(a,b,x,4097); //pow(2,mantissalength/2) + 1;
    }

    static double multiplication_error(double a, double b, double x) {
      return multiplication_error(a,b,x,67108865); //pow(2,mantissalength/2) + 1;
    }

    static T division_error(T a, T b, T x) {
      return (x*b - a - multiplication_error(x, b, x*b)) / b;
    }

  public:
    freal<T>(T value, T error) : value(value), error(error) { }
    freal<T>(T value) : value(value), error((T)0) { }
    freal<T>() : value((T)0), error((T)0) { }
    void operator+=(freal<T> rhs) {
      T value = this->value + rhs.value;
      T localerror = addition_error(this->value,rhs.value,value);
      this->value = value;
      this->error += rhs.error + localerror;
    }
    void operator-=(freal<T> rhs) {
      T value = this->value - rhs.value;
      T localerror = addition_error(this->value,-rhs.value,value);
      this->value = value;
      this->error -= rhs.error + localerror;
    }
    void operator*=(freal<T> rhs) {
      T value = this->value * rhs.value;
      T localerr = multiplication_error(this->value,rhs.value,value);
      this->error = this->value*rhs.error+this->error*rhs.value + localerr;
      this->value = value;
    }
    void operator/=(freal<T> rhs) {
      T recip = 1.0 / rhs.value;
      T newval = this->value * recip;
      T localerr = division_error_fl(this->value,rhs.value,newval);
      this->value = newval;
      this->error = recip*this->error - recip*newval*rhs.error + localerr;
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
  friend freal<T1> sin(const freal<T1> &);
  template <typename T1>
  friend freal<T1> cos(const freal<T1> &);
  template <typename T1>
  friend std::ostream& operator<<(std::ostream&, const freal<T1>&);
};

template<typename T>
freal<T> operator+(const freal<T> &g1,const freal<T> &g2){
  T value = g1.value + g2.value;
  T localerror = freal<T>::addition_error(g1.value,g2.value,value);
  T error = g1.error + g2.error + localerror;
  return freal<T>(value,error);
}

template<typename T>
freal<T> operator-(const freal<T> &g1,const freal<T> &g2){
  T value = g1.value - g2.value;
  T localerror = freal<T>::addition_error(g1.value,-g2.value,value);
  T error = g1.error - g2.error + localerror;
  return freal<T>(value,error);
}

template<typename T>
freal<T> operator-(const freal<T> &g1){
   return freal<T>(-g1.value,g1.error);
}

template<typename T>
freal<T> operator*(const freal<T> &g1,const freal<T> &g2){
   T value = g1.value * g2.value;
   T localerr = freal<T>::multiplication_error(g1.value,g2.value,value);
   T error = g1.value*g2.error+g1.error*g2.value + localerr;
   return freal<T>(value,error);
}

template<typename T>
freal<T> operator/(const freal<T> &g1,const freal<T> &g2){
   T recip,newval;
   recip = 1.0 / g2.value;
   newval = g1.value * recip;
   T localerr = freal<T>::division_error_fl(g1.value,g2.value,newval);
   return freal<T>(newval,recip*g1.error - recip*newval*g2.error + localerr);
}

template<typename T>
freal<T> sin(const freal<T> &g1){
  return afloat(sin(g1.value),cos(g1.value)*g1.error);
}

template<typename T>
freal<T> cos(const freal<T> &g1){
  return afloat(cos(g1.value),-sin(g1.value)*g1.error);
}

template<typename T>
std::ostream& operator<<(std::ostream &ost, const freal<T> &ad){
   ost << "[" << ad.value << ", " << ad.error << ", " << (ad.value - ad.error) << "]";
   return ost;
}

