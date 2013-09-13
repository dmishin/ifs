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

void Transform::set_scale( const point_t &s )
{
  t00 = s.x;
  t11 = s.y;
  t10 = t01 = 0;
}
point_t Transform::apply( const point_t& p)const
{
  return offset + 
    point_t( p.x*t00 + p.y*t01,
	     p.x*t10 + p.y*t11 );
}
void Transform::apply_inplace( point_t& p)const
{
  double xx = p.x*t00 + p.y*t01 + offset.x;
  double yy = p.x*t10 + p.y*t11 + offset.y;
  p.x = xx;
  p.y = yy;
}

Transform Transform::inverse()const
{
  double idet = 1.0/(t00*t11 - t10*t01);
  Transform invt;
  //Inverse 2x2 matrix
  invt.t00 = t11 * idet;
  invt.t11 = t00 * idet;
  invt.t01 = - t01 * idet;
  invt.t10 = - t10 * idet;

  //inverse offset
  invt.offset = point_t(0,0);
  invt.offset = invt.apply( -offset ); // o = T^-1 * (-o)

  return invt;
}

/**Multiply two transforms
   `this` is outer.
*/
Transform Transform::apply( const Transform&that )const
{
  // xt = o_this + T_this*(o_that + T_that * x)
  // xt = (o_this + T_this* o_that) + T_this * T_that * x;
  // xt = TFM_this(o_that) + 
  Transform mult;
  mult.offset = apply(that.offset);

  mult.t00 = t00 * that.t00 + t01*that.t10;
  mult.t01 = t00 * that.t01 + t01*that.t11;
  mult.t10 = t10 * that.t00 + t11*that.t10;
  mult.t11 = t10 * that.t01 + t11*that.t11;
  return mult;
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
  if (num_t < 0 || num_t > den_t) return false;
  
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
    if (segment_to_hrz_line_intersection( poly[i], 
					  poly[(i+1)%np], 
					  y, 
					  xs[i_intersect])){
      i_intersect ++;
    }
  }
  std::sort( xs, xs + i_intersect );
  double * p_end = std::unique( xs, xs + i_intersect );
  return p_end - xs;
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
