#include<cmath>

#include "../include/ran.h"

/*!
  \param[in] seed New seed
*/
void ran::setSeed(unsigned long long int seed) noexcept {
  v = 4101842887655102017LL;
  w = 1LL;
  u = seed^v; int64();
  v = u; int64();
  w = v; int64();
}

/*!
  \returns 64 bit random integer
*/
unsigned long long int ran::int64() noexcept {
  u = u * 2862933555777941757LL + 7046029254386353087LL;
  v ^= v >> 17; 
  v ^= v << 31; 
  v ^= v >> 8;
  w = 4294957665U*(w & 0xffffffff) + (w >> 32);
  unsigned long long int x = u^(u << 21); 
  x ^= x >> 35; 
  x ^= x << 4;
  return (x + v)^w;
}

/*!
  \param[in] minidx Minimum index
  \param[in] maxidx Maximum index
  \returns Random integer from the range [minidx,maxidx) -- that is,
     maxidx can't be generated
*/
unsigned int ran::selectFromRange(unsigned int minidx,
                                  unsigned int maxidx) noexcept {
  double relidx = (maxidx - minidx) * doub();
  return static_cast<unsigned int>(relidx) + minidx;
}

/*!
  \returns Gaussian deviate with mean 0 and variance 1
*/
double ran::gauss() noexcept {
  double u1,v1,x,y,q;
  do {
    u1 = doub();
    v1 = 1.7156 * (doub() - 0.5);
    x = u1 - 0.449871;
    y = fabs(v1) + 0.386595;
    q = x * x + y * (0.19600 * y - 0.25472 * x);
  } while (q > 0.27597 && (q > 0.27846 || v1 * v1 > -4. * log(u1) * u1 * u1));
  return v1 / u1;
}
