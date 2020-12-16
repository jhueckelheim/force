#ifndef __RANDOMHELPER
#define __RANDOMHELPER
template<typename T>
void print_binary(T number) {
  for(unsigned char i=0; i<sizeof(number)*8; i++) {
    putc((number & 1) ? '1' : '0', stdout);
    number = number >> 1;
  }
  putc('\n', stdout);
}

template<typename T>
T clearr(T in, char nbits) {
    return (in >> nbits) << nbits;
}

template<typename T>
T clearl(T in, char nbits) {
    return clearr(in, sizeof(in)*8-nbits)^in;
}

template<typename TI, int ML, int EL>
inline TI rand_fl() {
  static std::mt19937 mt{std::random_device{}()};
  static std::uniform_int_distribution<TI> mantissa_distribution;
  static std::uniform_int_distribution<int> exponent_distribution(-16,16);
  static std::uniform_int_distribution<TI> sign_distribution(0,1);
  TI mantissa = clearl<TI>(mantissa_distribution(mt), EL+1);
  TI exponent = (((TI)1<<EL-1)-1+exponent_distribution(mt))<<ML;
  TI sign = sign_distribution(mt)<<(ML+EL);
  print_binary<TI>(mantissa);
  print_binary<TI>(exponent);
  print_binary<TI>(sign);
  TI num = mantissa | exponent | sign;
  return num;
}

float rand_float() {
  unsigned int numi = rand_fl<unsigned int, 23, 8>();
  return *(reinterpret_cast<float*>(&numi));
}

double rand_double() {
  unsigned long long numi = rand_fl<unsigned long long, 52, 11>();
  return *(reinterpret_cast<double*>(&numi));
}

__float128 rand_quad() {
  unsigned long long numi[2];
  numi[0] = rand_fl<unsigned long long, 63, 0>();
  numi[1] = rand_fl<unsigned long long, 48, 15>();
  return *(reinterpret_cast<__float128*>(&(numi[0])));
}

long double rand_ldouble() {
  return (long double)rand_quad();
}
#endif
