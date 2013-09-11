#ifndef __IFS_INCLUDED__
#define __IFS_INCLUDED__
#include <vector>
#include "geometry.hpp"

typedef size_t pixel_index_t; //index of the pixel
typedef float amplitude_t; //amplitude of the pixel

struct source_pixel_ref_t{
  pixel_index_t i;
  amplitude_t k;
  source_pixel_ref_t(){};
  source_pixel_ref_t(pixel_index_t i_, amplitude_t k_):i(i_),k(k_){};
};

class MappingFunction{
public:
  virtual void map( const point_t *points, size_t n, point_t *out )=0;
};

class AffineMap: public MappingFunction{
public:
  point_t offset;
  double t00, t01, t10, t11;
  virtual void map( const point_t *points, size_t n, point_t *out );
};

class PixelMapping{
  friend class PixelMappingBuilder;
public:
  //Table of all pixel targets.
  typedef std::vector< source_pixel_ref_t > targets_array;
private:
  //Index of the first target for each pixel.
  //Last source is determined by the next index.
  //Contains one more value than number of pizels in the image
  std::vector< size_t > pixel_target_indices; 
  targets_array pixel_targets;
  //size of the image
  size_t width_, height_; 
protected:
  void add_target_link(size_t tgt, double k);
  void add_source_pixel(size_t src);
public:
  size_t width()const{ return width_; };
  size_t height()const{ return height_; };
  size_t n_relations()const{ return pixel_targets.size(); };
  PixelMapping( size_t w, size_t h );
  void get_targets_range(size_t idx, 
			 targets_array::const_iterator &ti, 
			 targets_array::const_iterator &ti_end)const;
};

class PixelMappingBuilder{
private:
  PixelMapping &mapping;
  MappingFunction *mapper;
  point_t top_left, scale, inv_scale;

  void register_pixel_image( size_t pixel_index, const std::vector<int> &pixels, size_t total_hits, int ix0, int iy0, int ix1, int iy1);

public:
  PixelMappingBuilder(PixelMapping &m);
  void setMapper( MappingFunction *m, const point_t &p0, const point_t &p1 );
  void build(size_t subpixels = 4);

};

/**Count, how many points inside rectangular box are inside the polygon.
   Box is given by 2 verties and n: number of points on the side.
   
 */
size_t count_points_inside_box( const point_t *poly, size_t np, const point_t &p0, const point_t &p1, size_t n);

class PixelMap;
void transform_pixel_map( const PixelMapping &mapping, const PixelMap &src, PixelMap &dst);

#endif
