#ifndef __PIXELMAP_HPP_INCLUDED__
#define __PIXELMAP_HPP_INCLUDED__
#include "pgm.hpp"
#include <vector>

class PixelMap{
public:
  typedef float pixel_t;
  std::vector<pixel_t> pixels;
  size_t width, height;

  PixelMap( size_t w, size_t h );
  void fill( pixel_t v );

  size_t pixel_idx(size_t x, size_t y)const{ return x+y*width; };
  pixel_t &pixel_ref(int x, int y){ return pixels[pixel_idx(x,y)]; };
  inline bool contains( int x, int y )const;
  pixel_t max_value()const;
  double norm2()const;
  void scale( double k );
  void normalize( pixel_t scale_to );
  void set_size(size_t w, size_t h){ width=w; height=h; pixels.resize(w*h);};
  void apply_gamma(double g);
  void apply_treshold(double h, double v);
  void swap(PixelMap &m);
  
};

bool PixelMap::contains( int x, int y )const
{
  return (x >= 0) && (y >=0) && (x < (int)width) && (y < (int)height);
}


class PixelMapWriter: public MonochromeImageWriter{
  PixelMap &pixmap;
  size_t cur_pixel_index;
public:
  PixelMapWriter( PixelMap &p):pixmap(p),cur_pixel_index(0){};
  virtual void set_size(size_t w, size_t h){ 
    pixmap.set_size(w,h); 
    cur_pixel_index = 0;
  };
  virtual void write_pixels(unsigned char *pixels, size_t n_pixels);
};

class PixelMapReader: public MonochromeImageReader{
  PixelMap &pixmap;
  size_t cur_pixel_index;
public:
  PixelMapReader(PixelMap &p):pixmap(p),cur_pixel_index(0){};
  virtual size_t read_pixels(unsigned char *pixels, size_t n_pixels);
  virtual size_t width(){ return pixmap.width; };
  virtual size_t height(){ return pixmap.height; };
};
#endif
