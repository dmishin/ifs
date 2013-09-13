#include<iostream>
#include<cassert>
#include<cmath>
#include<stdexcept>
#include<string>
#include <fstream>
#include <algorithm>

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


template<typename T1, typename T2>
struct static_caster{ T2 operator()(const T1 &x)const{ return static_cast<T2>(x); }; };
int main(int argc, char *argv[])
{
  using namespace std;
  try{
  if (false){ //test rendering polygon
    cout<<"Tesing render poly"<<endl;
    point_t poly[]={point_t(0,0), point_t(1,2),point_t(3,2.1), point_t(4,0.5)};
    vector<int> pix(100*100);
    render_polygon_aa_scanline(poly,4,
			       pix, 0,0, 100, 100, 15, 
			       point_t(4.0/100, 4.0/100), point_t(0,0));
    PixelMap img(100,100);
    transform(pix.begin(),pix.end(),img.pixels.begin(), static_caster<int,double>() );
    img.scale(1.0/255.0);
    PixelMapReader r(img);
    ofstream of("test-render.pgm",ios::binary);
    save_pgm(r, of);
    return 0;
  }
  


  //const char *in_file = "sample.pgm";
  const char *out_file = "test-tfm.pgm";
  //cout<<"Reading file "<<in_file<<endl;
  PixelMap in_image(800,800);
  in_image.fill(1.0);
  /*
  {
    PixelMapWriter w(in_image);
    ifstream in_fstream(in_file, ios::binary);
    read_pgm(in_fstream, w);
  }
  cout<<"Read image "<<in_file<<" :"<<in_image.width<<"x"<<in_image.height<<"px"<<endl;
  */
  //run simple testst
  test_is_to_the_right();
  test_segment_intersects_hrz_ray();
  test_is_inside_polygon();
  test_count_points_inside_box();

  //sample pixel mapping
  AffineMap map;
  map.tfm.rot_scale(-3.14159268/4, 0.5);
  map.tfm.offset = point_t(-0.2,0);

  AffineMap map1;
  map1.tfm.rot_scale(-3.14159268/8, 0.7);
  map1.tfm.offset = point_t(0.2,0);

  cout<<"Building mapping..."<<endl;
  cout.flush();
  PixelMapping mapping(in_image.width, in_image.height);
  {
    PixelMappingBuilder builder(mapping);
    builder.setMapper( &map, point_t(-1,-1), point_t(1,1) );
    builder.build( 4 );
  }
  cout<<"Mapping 1 built."
	   <<"Total relations:"<<mapping.n_relations()
	   <<endl;

  PixelMapping mapping1(in_image.width, in_image.height);
  {
    PixelMappingBuilder builder(mapping1);
    builder.setMapper( &map1, point_t(-1,-1), point_t(1,1) );
    builder.build( 4 );
  }
  cout<<"Mapping 2 built."
	   <<"Total relations:"<<mapping1.n_relations()
	   <<endl;

  size_t n_iters = 200;
  cout<<"Transforming for "<<n_iters<<" iterations ..."<<endl;
  PixelMap out_image(in_image.width,in_image.height);
  for( size_t iter=0; iter < n_iters; ++iter){
	  cout<<"."; cout.flush();
	  out_image.fill(0);
	  transform_pixel_map(mapping, in_image, out_image);
	  transform_pixel_map(mapping1, in_image, out_image);
	  out_image.normalize(1);
	  out_image.swap(in_image);
  }
  out_image.swap(in_image); //swap back
  out_image.apply_gamma(2);
  {
    PixelMapReader r(out_image);
    ofstream out_fstream(out_file, ios::binary);
    save_pgm(r, out_fstream);
  }
  cout<<"Saved "<<out_file<<endl;
  return 0;
  }catch(exception &e){
	  cerr<<"Exception raised:"<<e.what()<<endl;
	  return 1;
  }
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
