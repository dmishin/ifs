#include <cmath>
#include <stdexcept>
#include <algorithm>
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



/********************/

void bounds( point_t *pts, size_t n, point_t &a, point_t&b){
  a = b = pts[0];
  for(size_t i=1; i<n; ++i){
    a.x = std::min(a.x, pts[i].x);
    a.y = std::min(a.y, pts[i].y);
    b.x = std::max(b.x, pts[i].x);
    b.y = std::max(b.y, pts[i].y);
  }
}

/**True, if point x is to the right of the ray a->b*/
bool is_to_the_right(const point_t&a, const point_t&b, const point_t&x)
{
  point_t dv = b-a;
  point_t dp = x-a;
  return dp.x*dv.y - dp.y*dv.x >= 0;
}


bool segment_intersects_hrz_ray(const point_t&a, const point_t&b, const point_t&p0)
{
  // dx*t + x0 = k+px
  // dy*t + y0 = py
  // t in [0,1); k >= 0
  //
  // t = (py - y0)/dy
  // k = dx*t + (x0 - px)
  point_t d = b-a;
  if (d.y == 0.0) return false;
  //double t = (p0.y-a.y) / d.y;
  double num_t = p0.y-a.y, den_t = d.y;
  if (den_t < 0){ num_t = -num_t; den_t = -den_t; };
  if (num_t < 0 || num_t >= den_t) return false;
  //double k = a.x - p0.x + d.x*t;
  double num_k = (a.x - p0.x)*den_t +d.x*num_t;
  //double den_t = num_t;
  return num_k >= 0;
  //DONE: optimize without using division;
}


/**Intersection is inclusive
*/
bool segment_to_hrz_line_intersection(const point_t&a, const point_t&b, double y, double &x)
{
  point_t d = b-a;
  if (d.y == 0.0) return false; //horizontal segments has no intersections

  //double t = (p0.y-a.y) / d.y;
  double num_t = y-a.y, den_t = d.y;
  if (den_t < 0){ 
    num_t = -num_t; 
    den_t = -den_t; 
  };
  if (num_t <= 0 || num_t >= den_t) return false;
  
  x = d.x * num_t/den_t + a.x;
  return true;
}



/**builds list of intersection points between polygon and horizontal line.
   Points are sorted, if any.
   returns number of points found
 */
size_t polygon_to_hrz_line_intersections( const point_t* poly, size_t np, double y, double *xs, size_t max_xs)
{
  size_t i_intersect=0;
  for(size_t i=0; i<np && i_intersect < max_xs; ++i){
    if (segment_to_hrz_line_intersection( poly[i], poly[(i+1)%np], y, xs[i_intersect])){
      i_intersect ++;
    }
  }
  std::sort( xs, xs + i_intersect );
  return i_intersect;
}


bool is_inside_polygon(const point_t *poly, size_t n_points, const point_t &p)
{
  bool is_inside = false;
  for(size_t i=0; i<n_points; ++i){
    if (segment_intersects_hrz_ray(poly[i],poly[(i+1)%n_points],p))
      is_inside = ! is_inside;
  }
  return is_inside;
}
