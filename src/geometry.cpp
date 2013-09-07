#include <cmath>
#include <stdexcept>
#include "geometry.hpp"
#include "util.hpp"

double& Transform::as_vector( size_t idx )
{
  switch(idx){
  case 0: return offset.x;
  case 1: return offset.y;
  case 2: return t00;
  case 3: return t01;
  case 4: return t10;
  case 5: return t11;
  default: throw std::logic_error("bad index");
  }
}

double Transform::distance( const Transform &t )const
{
  double s=0;
  for(size_t i=0; i<6; ++i){
    s += sqr( as_vector(i)-t.as_vector(i));
  }
  return s;
}

void Transform::rot_scale( double alpha, double s)
{
  t00 = t11 = s*cos(alpha);
  double sn = sin(alpha)*s;
  t01 = sn;
  t10 = -sn;
}

point_t Transform::apply( const point_t& p)const
{
  return offset + 
    point_t( p.x*t00 + p.y*t01,
	     p.x*t10 + p.y*t11 );
}

