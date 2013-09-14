#ifndef __UTIL_HPP_INCLUDED__
#define __UTIL_HPP_INCLUDED__

#include <algorithm>

inline double sqr(double x){ return x*x; };
double random_double();

template<typename T>
T cap( T a, T b, T x)
{
  using namespace std;
  return min(max(x,a),b);
}

#endif
