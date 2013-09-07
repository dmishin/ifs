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

#endif
