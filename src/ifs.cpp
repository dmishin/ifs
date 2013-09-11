//reasonably fast IFS (iterfated functions system) fractal compuatation

#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include<iostream>
#include "geometry.hpp"
#include "ifs.hpp"
#include "util.hpp"
std::ostream &operator <<(std::ostream & os, const point_t &p)
{
  return os <<"["<<p.x<<","<<p.y<<"]";
}

void AffineMap::map( const point_t *points, size_t n, point_t *out )
{
  for(size_t i=0; i<n; ++i){
    out[i] = offset + point_t(t00*points[i].x + t01*points[i].y,
			      t10*points[i].x + t11*points[i].y);
  }
}


PixelMapping::PixelMapping(size_t w, size_t h)
  :width_(w), height_(h)
{
  pixel_target_indices.resize(width_*height_);
}

void PixelMapping::get_targets_range(
  size_t idx, 
  PixelMapping::targets_array::const_iterator &ti,
  PixelMapping::targets_array::const_iterator &ti_end)
{
  ti = pixel_targets.begin() + pixel_target_indices[idx];
  if (idx + 1 >= pixel_target_indices.size())
    ti_end = pixel_targets.end();
  else
    ti_end = pixel_targets.begin() + pixel_target_indices[idx+1];
}

void PixelMappingBuilder::setMapper( MappingFunction *m, const point_t &p0, const point_t &p1 )
{
  mapper = m;
  top_left = p0;
  scale = (p1-p0)/point_t( mapping.width(), mapping.height());
  inv_scale = point_t( mapping.width(), mapping.height()) / (p1-p0);
}



/**Count, how many points inside rectangular box are inside the polygon.
   Box is given by 2 verties and n: number of points on the side.
   
 */
size_t count_points_inside_box( const point_t *poly, size_t np, const point_t &p0, const point_t &p1, size_t n)
{
  size_t inside = 0;
  point_t dp = (p1 - p0)*(1.0/n);
  point_t p;
  for(int iy=0; iy<(int)n; iy++){
    p.y = p0.y+dp.y*iy + dp.y*0.5;
    int left_bnd=n, right_bnd=n-1;
    for(int ix=0; ix<(int)n; ix++){
      p.x = p0.x + dp.x*ix + dp.x*0.5;
      if (is_inside_polygon( poly, np, p )){
	left_bnd = right_bnd = ix;
	break;
      }
    }
    for(int ix=n-1; ix>left_bnd;--ix){
      p.x = p0.x + dp.x*ix + dp.x*0.5;
      if (is_inside_polygon( poly, np, p )){
	right_bnd = ix;
	break;
      }
    }
    /*
    for(int ix=left_bnd; ix<=right_bnd;++ix){
      p.x = p0.x + dp.x*ix + dp.x*0.5;
      assert( is_inside_polygon( poly, np, p ));
    }
    */
    inside += (right_bnd - left_bnd + 1);
  }
  return inside;
}


void PixelMapping::add_source_pixel(size_t src)
{
  assert( src < pixel_target_indices.size() );
  pixel_target_indices[src] = pixel_targets.size();
}

void PixelMapping::add_target_link(size_t tgt, double k)
{
  pixel_targets.push_back(source_pixel_ref_t(tgt,k));
}

/**Anti-aliased polygon rendering
 * Result is array of intensities @poly
 * return value is total intensity (sum of array elements).
 */
size_t render_polygon_aa(const point_t *poly, size_t poly_points,
			 std::vector<int> &pixels, int ix0, int iy0, int ix1, int iy1, size_t subpixels, 
			 const point_t scale, const point_t top_left)
{
  int dest_w = ix1-ix0;
  int dest_h = iy1-iy0;
  pixels.resize(dest_w*dest_h); //row-by-row

  int total_hits = 0;
  if ((dest_w == 1) && (dest_h==1)){
    //special case: only one target pixel
    total_hits = 1;
    pixels[0] = 1;
  }else{
    size_t idx = 0;
    for( int yy=iy0; yy< iy1; ++yy){
      for( int xx=ix0; xx< ix1; ++xx, ++idx){
	point_t pp0 = point_t(xx,yy)*scale + top_left;
	point_t pp1 = pp0+scale;
	size_t hits = count_points_inside_box(poly,poly_points, pp0, pp1, subpixels);
	total_hits += hits;
	pixels[idx] = hits;
      }
    }
  }
  return total_hits;
}

size_t paint_segment_to_pixel_row( double x0, double x1, int *pixels, int size );
/**Anti-aliased polygon rendering
 * Result is array of intensities @poly
 * return value is total intensity (sum of array elements).
 * scanline-based version
 */
size_t render_polygon_aa_scanline(const point_t *poly, size_t poly_points,
			 std::vector<int> &pixels, int ix0, int iy0, int ix1, int iy1, size_t subpixels, 
			 const point_t scale, const point_t top_left)
{
  int dest_w = ix1-ix0;
  int dest_h = iy1-iy0;
  pixels.resize(dest_w*dest_h); //row-by-row

  const size_t max_intersections = 4;//4 is the maximum possible value, for non-convex images (rare)
  double intersections[max_intersections];

  int total_hits = 0;
  if ((dest_w == 1) && (dest_h==1)){
    //special case: only one target pixel
    total_hits = 1;
    pixels[0] = 1;
  }else{
    std::fill(pixels.begin(), pixels.end(), 0);

    for( int yy=iy0; yy< iy1; ++yy){
      double y_base = top_left.y + scale.y*yy;
      for ( size_t iy = 0; iy < subpixels; ++iy ){
	double y = y_base + (iy+0.5) * scale.y / subpixels;
	size_t n_inters = polygon_to_hrz_line_intersections( poly, poly_points, y, intersections, max_intersections );
	if (n_inters % 2 == 1){
	  std::cerr<<"warning: odd number of intersections"<<std::endl;
	  n_inters --; //this should not be... but just in case - ensure that we have 
	}
	//now for each segment, "paint" it to the pixelmap
	for( size_t j=0; j<n_inters; j += 2 ){
	  total_hits += paint_segment_to_pixel_row( (intersections[j]-top_left.x)/scale.x-ix0, 
						    (intersections[j+1]-top_left.x)/scale.x-ix0, 
						    &pixels[(yy-iy0)*dest_w], dest_w);
	}
      }
    }
  }
  return total_hits;
}

size_t paint_segment_to_pixel_row( double x0, double x1, int *pixels, int size )
{
  assert( x1 >= x0 );
  if (x1 < 0) return 0;
  if (x0 > size) return 0;

  double cx0 = cap( 0.0, (double)size, x0);
  double cx1 = cap( 0.0, (double)size, x1);
  int i0 = (int)floor(cx0);
  int i1 = (int)ceil(cx1);

  assert( i0 >= 0);
  assert( i0 <= size);  
  assert( i1 >= 0);  
  assert( i1 <= size);

  double p0 = cx0 - floor(cx0);
  double p1 = ceil(cx1) - cx1;

  for (int i=i0; i<i1; ++i){
    pixels[i] += 256;
  }
  //first and last pixels are incomplete: decrease them.
  int d0 = (int)(256 * p0);
  int d1 = (int)(256 * p1);
  pixels[i0] -= d0;
  pixels[i1-1] -= d1;

  return 256 * (i1-i0) - d0 - d1; //total sum of added values
}

void fill_hrz_row_of_points( point_t *points, size_t npoints, const point_t &p0, double dx)
{
  for(size_t i=0; i<npoints;++i){
    points[i].y = p0.y;
    points[i].x = p0.x + dx * i;
  }
}

void PixelMappingBuilder::build(size_t subpixels)
{
  size_t width = mapping.width(), height = mapping.height();
  assert(mapper);
  std::vector<int> dest_pixmap;
  point_t dst[4];
  point_t p0, p1;
  std::vector< point_t > row_src(width+1);
  std::vector< point_t > cur_row_dst(width+1), next_row_dst(width+1);
  
  for(int y=0; y<(int)height; ++y){
    //Transform images of the pixels. Do it row-by-row
    if (y == 0){
      fill_hrz_row_of_points( &row_src[0], row_src.size(), top_left+point_t(0,y)*scale, scale.x );
      mapper->map(&row_src[0], row_src.size(), &cur_row_dst[0]);
    }else{ //reuse value from the previous step
      cur_row_dst.swap(next_row_dst);
    }
    fill_hrz_row_of_points( &row_src[0], row_src.size(), top_left+point_t(0,y+1)*scale, scale.x );
    mapper->map(&row_src[0], row_src.size(), &next_row_dst[0]);

    //Now, for each pixel in the row, build its anti-aliased image and store it in the mapping
    for(int x=0; x<(int)width; ++x){
      size_t pixel_index = (size_t)(x + y*width);
      //get list of destination pixels
      dst[0] = cur_row_dst[x];
      dst[1] = next_row_dst[x];
      dst[2] = next_row_dst[x+1];
      dst[3] = cur_row_dst[x+1];

      bounds( dst, 4, p0, p1);
      point_t scr_p0 = (p0 - top_left) * inv_scale;
      point_t scr_p1 = (p1 - top_left) * inv_scale;

      int ix0 = cap(0, (int)width-1, (int)floor(scr_p0.x));
      int iy0 = cap(0, (int)height-1, (int)floor(scr_p0.y));
      int ix1 = cap(0, (int)width-1, (int)ceil(scr_p1.x));
      int iy1 = cap(0, (int)height-1, (int)ceil(scr_p1.y));
      int total_hits = render_polygon_aa_scanline(dst, 4, dest_pixmap, ix0, iy0, ix1, iy1, subpixels, scale, top_left);
      //Now we have an anti-aliased map of target pixels for the given source pixel.
      //Record the source-target relation
      register_pixel_image( pixel_index, dest_pixmap, total_hits, ix0, iy0, ix1, iy1);
    }
  }
}

/**Take grayscale antialiased image of a transformed pixel, and add it to the mapping
 */
void PixelMappingBuilder::register_pixel_image( size_t pixel_index, const std::vector<int> &pixels, size_t total_hits, int ix0, int iy0, int ix1, int iy1)
{
  if( total_hits == 0) return;
  double kpix = 1.0/double(total_hits);

  mapping.add_source_pixel(pixel_index);

  size_t idx=0;
  for( int yy=iy0; yy< iy1; ++yy){
    for( int xx=ix0; xx< ix1; ++xx, ++idx){
      if (pixels[idx] > 0){
	double k = double(pixels[idx])*kpix;
	mapping.add_target_link( yy*mapping.width()+xx, k);
      }
    }
  }
}



PixelMappingBuilder::PixelMappingBuilder(PixelMapping &m)
  :mapping(m)
  ,mapper(NULL)
{
}
