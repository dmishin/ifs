#include<iostream>
#include<cassert>
#include<cmath>
#include "ifs.hpp"

void test_segment_intersects_hrz_ray()
{
  assert( segment_intersects_hrz_ray(point_t(0,0), point_t(1,1), 
				     point_t(0,0.5)));
  assert( segment_intersects_hrz_ray(point_t(1,1), point_t(0,0), 
				     point_t(0,0.5)));
  assert( segment_intersects_hrz_ray(point_t(0,0), point_t(1,1), 
				     point_t(-1,0.5)));
  assert(!segment_intersects_hrz_ray(point_t(0,0), point_t(1,1), 
				     point_t(1,0.5)));
}

void test_count_points_inside_box()
{
  point_t poly[4]={point_t(0,0), point_t(0,1), point_t(1,1),point_t(1,0)};
  size_t np = 4;
  assert( count_points_inside_box(poly, np, point_t(0.1,0.1), point_t(0.9, 0.9), 5) == 25 );

  assert( count_points_inside_box(poly, np, point_t(0.1,0.1), point_t(1.9, 0.9), 5) < 25 );
  assert( count_points_inside_box(poly, np, point_t(0.1,0.1), point_t(1.9, 0.9), 5) > 0 );

  assert( count_points_inside_box(poly, np, point_t(1.1,0.1), point_t(1.9, 0.9), 5) == 0 );
  
}


void test_is_to_the_right()
{
  assert(is_to_the_right( point_t(0,0), point_t(0,1), point_t(0.5,0.5)));
  assert(is_to_the_right( point_t(0,0), point_t(0,1), point_t(1,0)));
  assert(is_to_the_right( point_t(0,0), point_t(0,1), point_t(0.5,-0.5)));
  assert(!is_to_the_right( point_t(0,0), point_t(0,1), point_t(-0.5,0.5)));
  assert(!is_to_the_right( point_t(0,0), point_t(0,1), point_t(-1,0)));
}

void test_is_inside_polygon()
{
  point_t poly[4]={point_t(0,0), point_t(0,1), point_t(1,1),point_t(1,0)};
  point_t ipoly[4]={point_t(0,0), point_t(0,1),point_t(1,0), point_t(1,1)};
  size_t np = 4;
  assert(is_inside_polygon(poly,np, point_t(0.5,0.5)));
  assert(is_inside_polygon(ipoly,np, point_t(0.5,0.5)));
  assert(is_inside_polygon(poly,np, point_t(0.1,0.5)));
  assert(is_inside_polygon(poly,np, point_t(0.5,0.7)));
  assert(!is_inside_polygon(poly,np, point_t(1.5,0.5)));
  assert(!is_inside_polygon(poly,np, point_t(-.5,0.5)));
  assert(!is_inside_polygon(poly,np, point_t(.5,1.5)));
  assert(!is_inside_polygon(poly,np, point_t(1.5,-0.5)));
  assert(!is_inside_polygon(poly,np, point_t(1.5,1.5)));
}

int main(int argc, char *argv[])
{
  //run simple testst
  test_is_to_the_right();
  test_segment_intersects_hrz_ray();
  test_is_inside_polygon();
  test_count_points_inside_box();

  //sample pixel mapping
  AffineMap map;
  map.offset = point_t(0,0);
  double alpha = 0.4;
  double s = 3;
  map.t11 = map.t00 = cos(alpha)*s;
  map.t01 = sin(alpha)*s;
  map.t10 = -sin(alpha)*s;

  std::cout<<"Building mapping..."<<std::endl;
  std::cout.flush();
  PixelMapping mapping(2000, 2000);
  mapping.setMapper( &map, point_t(-1,-1), point_t(1,1) );
  mapping.build( 2 );
  std::cout<<"Mapping built."
	   <<"Total relations:"<<mapping.pixel_targets.size()
	   <<std::endl;
  std::cout<<"Targets for the point 150,100:"<<std::endl;
  size_t idx = 1050 + 1000 * mapping.width;
  size_t ti, ti_end;
  mapping.get_targets_range(idx, ti, ti_end);
  double sk = 0;
  for(size_t i=ti; i!= ti_end; ++i){
    size_t j = mapping.pixel_targets[i].i;
    std::cout<<"("<<j%mapping.width<<","<<j/mapping.width<<") -- "<<mapping.pixel_targets[i].k
	     <<std::endl;
    sk += mapping.pixel_targets[i].k;
  }
  std::cout<<"  total k: "<<sk<<std::endl;
  return 0;
}
