#include <iostream>
#include <math.h>

template<typename T>
class ForceFloat {
  private:
    T value, error;
    T addition_error(T a, T b, T x) {
        if (fabs(a) >= fabs(b)) {
            return (x - a) - b;
        } else {
            return (x - b) - a;
        }
    }

    T multiplication_error(T a, T b, T x, int m) {
      T au, al, bu, bl;
      au = (a - a*m) + a*m;
      al = a - au;
      bu = (b - b*m) + b*m;
      bl = b - bu;
      return x - au*bu - (au*bl + al*bu) - al*bl;
    }

    float multiplication_error(float a, float b, float x) {
      return multiplication_error(a,b,x,4097); //pow(2,mantissalength/2) + 1;
    }

    double multiplication_error(double a, double b, double x) {
      return multiplication_error(a,b,x,67108865); //pow(2,mantissalength/2) + 1;
    }

  public:
    ForceFloat<T>(T value, T error) : value(value), error(error) { }
    ForceFloat<T>(T value) : value(value), error((T)0) { }

    ForceFloat<T> operator+(const ForceFloat<T> &g1,const ForceFloat<T> &g2){
      T value = g1.value + g2.value;
      T localerror = addition_error(g1.value,g2.value,value);
      T error = g1.error + g2.error + localerror;
      return ForceFloat<T>(value,error);
    }

    ForceFloat<T> operator-(const ForceFloat<T> &g1,const ForceFloat<T> &g2){
      T value = g1.value - g2.value;
      T localerror = addition_error(g1.value,-g2.value,value);
      T error = g1.error - g2.error + localerror;
      return ForceFloat<T>(value,error);
    }

    ForceFloat<T> operator-(const ForceFloat<T> &g1){
       return ForceFloat<T>(-g1.value,g1.error);
    }

    ForceFloat<T> operator*(const ForceFloat<T> &g1,const ForceFloat<T> &g2){
       T value = g1.value * g2.value;
       T localerr = multiplication_error_fl(g1.value,g2.value,value);
       T error = g1.value*g2.error+g1.error*g2.value + localerr;
       return ForceFloat<T>(value,error);
    }

    ForceFloat<T> operator/(const ForceFloat<T> &g1,const ForceFloat<T> &g2){
       T recip,newval;
       recip = 1.0 / g2.value;
       newval = g1.value * recip;
       T localerr = division_error_fl(g1.value,g2.value,newval);
       return ForceFloat<T>(newval,recip*g1.error - recip*newval*g2.error + localerr);
    }

  friend std::ostream& operator<<(std::ostream&, const ForceFloat<T>&);
  friend ForceFloat<T> sin(const ForceFloat<T> &);
  friend ForceFloat<T> cos(const ForceFloat<T> &);
}

template<typename T>
ForceFloat<T> sin(const ForceFloat<T> &g1){
  return afloat(sin(g1.value),cos(g1.value)*g1.error);
}

template<typename T>
ForceFloat<T> cos(const ForceFloat<T> &g1){
  return afloat(cos(g1.value),-sin(g1.value)*g1.error);
}

template<typename T>
std::ostream& operator<<(std::ostream &ost, const ForceFloat<T> &ad){
   ost << "[" << ad.value << ", " << ad.error << ", " << (ad.value - ad.error) << "]";
   return ost;
}

