#ifndef __GEOMETRY_HPP_INCLUDED__
#define __GEOMETRY_HPP_INCLUDED__

struct point_t{
  double x;
  double y;
  point_t(){};
  point_t(double x_, double y_):x(x_), y(y_){};
  point_t operator +(const point_t& p)const{
    return point_t(x+p.x, y+p.y);
  };
  point_t operator -(const point_t& p)const{
    return point_t(x-p.x, y-p.y);
  };
  point_t operator *(const point_t& k)const{
    return point_t(x*k.x, y*k.y);
  };
  point_t operator *(double k)const{
    return point_t(x*k, y*k);
  };
  point_t &operator +=(const point_t &p){
    x+=p.x; y+=p.y;
    return *this;
  };
};


class Transform{
public:
  double t00,t01, t10, t11;
  point_t offset;
  point_t apply( const point_t& p)const;
  void rot_scale( double alpha, double s);
  double distance( const Transform &t )const;
  double &as_vector( size_t idx );
  double as_vector( size_t idx)const{ 
    return const_cast<Transform*>(this)->as_vector(idx);
  }
};

void bounds( point_t *pts, size_t n, point_t &a, point_t&b);
/**True, if point x is to the right of the ray a->b*/
bool is_to_the_right(const point_t&a, const point_t&b, const point_t&x);
bool segment_intersects_hrz_ray(const point_t&a, const point_t&b, const point_t&p0);
/**Intersection is inclusive
*/
bool segment_to_hrz_line_intersection(const point_t&a, const point_t&b, double y, double &x);
/**Intersection is inclusive
*/
bool segment_to_hrz_line_intersection(const point_t&a, const point_t&b, double y, double &x);
/**builds list of intersection points between polygon and horizontal line.
   Points are sorted, if any.
   returns number of points found
 */
size_t polygon_to_hrz_line_intersections( const point_t* poly, size_t np, double y, double *xs, size_t max_xs);
bool is_inside_polygon(const point_t *poly, size_t n_points, const point_t &p);

#endif
