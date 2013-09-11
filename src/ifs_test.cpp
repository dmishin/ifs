#include<iostream>
#include<cassert>
#include<cmath>
#include<stdexcept>
#include<string>
#include <fstream>

#include "ifs.hpp"
#include "pixelmap.hpp"
#include "pgm.hpp"

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
  const char *in_file = "test.pgm";
  const char *out_file = "test-tfm.pgm";
  PixelMap in_image(0,0);
  {
    PixelMapWriter w(in_image);
    std::ifstream in_fstream(in_file, std::ios::binary);
    read_pgm(in_fstream, w);
  }
  std::cout<<"Read image "<<in_image.width<<"x"<<in_image.height<<std::endl;
  
  //run simple testst
  test_is_to_the_right();
  test_segment_intersects_hrz_ray();
  test_is_inside_polygon();
  test_count_points_inside_box();

  //sample pixel mapping
  AffineMap map;
  map.offset = point_t(0,0);
  double alpha = 0.3;
  double s = 2;
  map.t11 = map.t00 = cos(alpha)*s;
  map.t01 = sin(alpha)*s;
  map.t10 = -sin(alpha)*s;

  std::cout<<"Building mapping..."<<std::endl;
  std::cout.flush();
  PixelMapping mapping(in_image.width, in_image.height);
  {
    PixelMappingBuilder builder(mapping);
    builder.setMapper( &map, point_t(-1,-1), point_t(1,1) );
    builder.build( 8 );
  }
  std::cout<<"Mapping built."
	   <<"Total relations:"<<mapping.n_relations()
	   <<std::endl;
  std::cout<<"Transforming..."<<std::endl;
  PixelMap out_image(0,0);
  transform_pixel_map(mapping, in_image, out_image);
  out_image.normalize(1);
  {
    PixelMapReader r(out_image);
    std::ofstream out_fstream(out_file, std::ios::binary);
    save_pgm(r, out_fstream);
  }
  std::cout<<"Saved "<<out_file<<std::endl;
  return 0;
}



  /*
  std::cout<<"Targets for the point 150,100:"<<std::endl;
  size_t idx = 1050 + 1000 * mapping.width();
  PixelMapping::targets_array::const_iterator ti, ti_end;
  mapping.get_targets_range(idx, ti, ti_end);
  double sk = 0;
  for(PixelMapping::targets_array::const_iterator i=ti; i!= ti_end; ++i){
    size_t j = i->i;
    size_t w = mapping.width();
    std::cout<<"("<<j%w<<","<<j/w<<") -- "<<i->k
	     <<std::endl;
    sk += i->k;
  }
  std::cout<<"  total k: "<<sk<<std::endl;
  */
