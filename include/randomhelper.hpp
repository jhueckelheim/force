#ifndef __RANDOMHELPER
#define __RANDOMHELPER
#include <random>

/*
 * A slightly unusual random number generator for floating point numbers:
 * Standard generators may not support the full long-double or quad
 * resolution, meaning that they generate numbers only in lower precision
 * and leave some mantissa bits non-random. Also, they tend to be defined
 * as a distribution in real numbers, which makes some floating point
 * numbers more likely to be chosen than others during rounding/truncation,
 * due to their uneven spacing.
 * To fix this, we build floating point numbers from scratch, using a
 * random bit sequence as the mantissa, a random bit for the sign, and a
 * random number from the allowable exponent range for the exponent.
 * The sign|exponent|mantissa bit sequences are then combined and cast
 * into a floating point number.
 * */

// Helper function to print the bit sequence of a given value to stdout.
// Used only during debugging.
template<typename T>
void print_binary(T number) {
  for(unsigned char i=0; i<sizeof(number)*8; i++) {
    putc((number & 1) ? '1' : '0', stdout);
    number = number >> 1;
  }
  putc('\n', stdout);
}

// Clear the upper-most `nbits` bits in `in`.
template<typename T>
T clearr(T in, char nbits) {
    return (in >> nbits) << nbits;
}

// Clear the lower-most `nbits` bits in `in`.
template<typename T>
T clearl(T in, char nbits) {
    return clearr(in, sizeof(in)*8-nbits)^in;
}

// Build a number of size `sizeof(TI)`, containing `ML` mantissa bits,
// `EL` exponent bits, and a sign bit. This function is used internally,
// there are user-friendly wrappers below for some standard number types.
template<typename TI, int ML, int EL>
inline TI rand_fl() {
  //static std::mt19937 mt{std::random_device{}()};
  static std::mt19937 mt;
  static std::uniform_int_distribution<TI> mantissa_distribution;
  static std::uniform_int_distribution<int> exponent_distribution(-16,16);
  static std::uniform_int_distribution<TI> sign_distribution(0,1);
  TI mantissa = clearl<TI>(mantissa_distribution(mt), EL+1);
  TI exponent = (((TI)1<<EL-1)-1+exponent_distribution(mt))<<ML;
  TI sign = sign_distribution(mt)<<(ML+EL);
  TI num = mantissa | exponent | sign;
  return num;
}

// Build a IEEE754 single precision number: 23 mantissa bits, 8 exponent
// bits, 1 sign bit, total size of an unsigned-int (32 bits).
float rand_float() {
  unsigned int numi = rand_fl<unsigned int, 23, 8>();
  return *(reinterpret_cast<float*>(&numi));
}

// Build a IEEE754 double precision number: 52 mantissa bits, 11 exponent
// bits, 1 sign bit, total size of an unsigned-long-long (64 bits).
double rand_double() {
  unsigned long long numi = rand_fl<unsigned long long, 52, 11>();
  return *(reinterpret_cast<double*>(&numi));
}

// Build a IEEE754 quad precision number: 112 mantissa bits, 15 exponent
// bits, 1 sign bit, total size of two unsigned-long-long (2 * 64 bits).
// Since there is no built-in integer type of 128bit length, we split
// the construction into two phases applied to consecutive 64bit locations.
// This means constructing a 48bit mantissa / 15bit exponent / 1bit sign
// and another 63bit mantissa / 0bit exponent / 1bit sign (which is then
// re-interpreted as a 64bit mantissa). Both stored next to each other
// and re-interpreted as a __float128 yields our quad precision number.
__float128 rand_quad() {
  unsigned long long numi[2];
  numi[0] = rand_fl<unsigned long long, 63, 0>();
  numi[1] = rand_fl<unsigned long long, 48, 15>();
  return *(reinterpret_cast<__float128*>(&(numi[0])));
}

// Build a long-double. Since this is not a IEEE754 number, things are a
// little weirder here. On many platforms, this is an 80bit number with
// no implicit upper mantissa bit, which complicates generation slightly.
// We instead just generate a __float128 and cast it to a long double.
// Down-casting should still yield a fully random long-double mantissa.
long double rand_ldouble() {
  return (long double)rand_quad();
}
#endif
