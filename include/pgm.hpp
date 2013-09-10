#ifndef __PGM__HPP_INCLUDED__
#define __PGM__HPP_INCLUDED__
#include <iostream>


struct MonochromeImageWriter{
  virtual void set_size(size_t w, size_t h) = 0;
  virtual void write_pixels(unsigned char *pixels, size_t n_pixels)=0;
};
struct MonochromeImageReader{
  virtual size_t width() = 0;
  virtual size_t height() = 0;
  virtual size_t read_pixels(unsigned char *pixels, size_t n_pixels)=0;
};

void save_pgm( MonochromeImageReader &pixels, std::ostream &ofile );
void read_pgm( std::istream &ifile, MonochromeImageWriter &pixels);

#endif
