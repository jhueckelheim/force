#include <iostream>
#include <math.h>

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
    freal<T>(T value, T error) : val(value), err(error) { }
    freal<T>(T value) : val(value), err((T)0) { }
    freal<T>() : val((T)0), err((T)0) { }
    T value() const { return this->val; }
    T error() const { return this->err; }
    T corrected_value() const { return this->val-this->err; }

    void operator+=(freal<T> rhs) {
      T value = this->val + rhs.val;
      T localerror = addition_error(this->val,rhs.val,value);
      this->val = value;
      this->err += rhs.err + localerror;
    }
    void operator-=(freal<T> rhs) {
      T value = this->val - rhs.val;
      T localerror = addition_error(this->val,-rhs.val,value);
      this->val = value;
      this->err -= rhs.err + localerror;
    }
    void operator*=(freal<T> rhs) {
      T value = this->val * rhs.val;
      T localerr = multiplication_error(this->val,rhs.val,value);
      this->err = this->val*rhs.err+this->err*rhs.val + localerr;
      this->val = value;
    }
    void operator/=(freal<T> rhs) {
      T recip = 1.0 / rhs.val;
      T newval = this->val * recip;
      T localerr = division_error(this->val,rhs.val,newval);
      this->val = newval;
      this->err = recip*this->err - recip*newval*rhs.err + localerr;
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
   recip = 1.0 / g2.val;
   newval = g1.val * recip;
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
std::ostream& operator<<(std::ostream &ost, const freal<T> &ad){
   ost << "[" << ad.value() << ", " << ad.error() << ", " << (ad.corrected_value()) << "]";
   return ost;
}

