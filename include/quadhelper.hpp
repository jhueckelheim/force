#ifndef _QUADHELPER
#define _QUADHELPER
typedef __float128 usereal;
usereal sqrt(usereal val) {
   return sqrtq(val);
}
std::ostream& operator<<(std::ostream &ost, const __float128 &ad){
   char buf[128];
   quadmath_snprintf(buf, sizeof(buf), "%.36Qg", ad);
   ost<<buf;
   return ost;
}
#endif
