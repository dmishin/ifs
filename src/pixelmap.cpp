#include <cmath>
#include "pixelmap.hpp"

void PixelMap::normalize(PixelMap::pixel_t scale_to)
{
  scale(double(scale_to)/double(max_value()));
}

void PixelMap::apply_treshold(double h, double v)
{
	std::vector<pixel_t>::iterator i, e;
	for(i=pixels.begin(), e=pixels.end(); i != e; ++i ){
		*i = (*i >= h)? v : 0;
	}
}

void PixelMap::apply_gamma( double g )
{
  double p = 1.0/g;
  for(size_t i=0; i<pixels.size(); ++i){
    pixels[i]=(pixel_t)pow((double)pixels[i], p);
  }  
}
void PixelMap::scale( double k )
{
  for(size_t i=0; i<pixels.size(); ++i){
    pixels[i]=(pixel_t)(pixels[i]*k);
  }
}
PixelMap::PixelMap( size_t w, size_t h )
{
  width=w;
  height=h;
  pixels.resize(w*h);
  fill(0);
}
void PixelMap::fill( pixel_t v )
{
  for(size_t i=0; i<pixels.size(); ++i){
    pixels[i]=v;
  }
}

PixelMap::pixel_t PixelMap::max_value()const
{
  pixel_t v=0;
  for(size_t i=0; i<pixels.size(); ++i){
    if (v<pixels[i]) v=pixels[i];
  }
  return v;
}


void PixelMapWriter::write_pixels(unsigned char *pixels, size_t n_pixels)
{
  for(size_t i=0; i<n_pixels; ++i){
    pixmap.pixels[cur_pixel_index + i] = 
      PixelMap::pixel_t((double)pixels[i] / (double)255.0);
  }
  cur_pixel_index += n_pixels;
}

size_t PixelMapReader::read_pixels(unsigned char *pixels, size_t n_pixels)
{
  size_t n_left = pixmap.pixels.size() - cur_pixel_index;
  size_t n_read = std::min(n_pixels, n_left);
  for(size_t i=0; i<n_read; ++i){
    pixels[i] = 
      (unsigned char)(int)(255 * pixmap.pixels[cur_pixel_index + i]);
  }
  cur_pixel_index += n_read;
  return n_read;
}


void PixelMap::swap(PixelMap &m)
{
	m.pixels.swap(pixels);
	std::swap(width, m.width);
	std::swap(height, m.height);
}
