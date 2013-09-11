#ifndef __UTIL_HPP_INCLUDED__
#define __UTIL_HPP_INCLUDED__

#include <algorithm>

inline double sqr(double x){ return x*x; };
double random_double();

template<typename T>
T cap( const T&a, const T&b, const T&x)
{
  return std::min(std::max(x,a),b);
}


#endif
