#include "pgm.hpp"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include "ifs-rand.hpp"

size_t read_pgm_header_int(std::istream &ifile);
void skip_to_eol(std::istream &ifile);


void save_pgm( const PixelMap &pixels, std::ostream &ofile, double gamma/*=1.0*/ )
{
  ofile << "P5\n"; //binary graymap
  ofile << pixels.width << " " << pixels.height << "\n";
  ofile << "255\n";
  PixelMap::pixel_t max_value = pixels.max_value();
  if (max_value < 1) max_value=1;
  double igamma = 1.0/gamma;
  for( size_t i=0; i<pixels.pixels.size(); ++i){
    double v = (double)pixels.pixels[i]/(double)max_value;
    v = pow(v, igamma);
    unsigned char g = (unsigned char)(int)(v*255);
    ofile << (char)g;
  };
}

void read_pgm_header( std::istream &ifile, size_t &w, size_t &h, size_t &maxcolor )
{
  char magic[2];
  ifile.read(&magic[0], 2);
  if (magic[0] !='P' || magic[1] != '5')
    throw std::logic_error("Bad magic number");
  w = read_pgm_header_int(ifile);
  h = read_pgm_header_int(ifile);
  maxcolor = read_pgm_header_int(ifile);
}

void skip_to_eol(std::istream &ifile)
{
  char c;
  while(true){
    ifile.read( &c, 1);
    if (c=='\n' || c=='\r') return;
  }
}
size_t read_pgm_header_int(std::istream &ifile)
{
  char c;
  size_t x = 0;
  bool has_number = false;
  while(true){
    ifile.read(&c,1);
    if (c==' ' || c=='\t' || c=='\r' || c=='\n'){
      if (has_number)
	break;
      else
	continue;
    }
    if (c >= '0' && c <= '9'){
      has_number=true;
      x = x*10 + (size_t)(c - '0');
      continue;
    }
    if (c == '#'){
      skip_to_eol(ifile);
      continue;
    }
    throw std::logic_error( "unexpected character in PGM file");
  }
  return x;
}
PixelMap * read_pgm( std::istream &ifile )
{
  size_t w, h, maxcolor;
  read_pgm_header( ifile, w, h, maxcolor );
  if (maxcolor > 255)
    throw std::logic_error( "only 8-bit supported");
  if (w ==0 || h==0 )
    throw std::logic_error( "zero size of image");
  std::cerr<<"Reading image "<<w<<"x"<<h<<" with "<<maxcolor<<" colors\n";
  PixelMap *pix = new PixelMap(w,h);
  for(size_t y=0; y<h; ++y){
    for(size_t x=0; x<w; ++x){
      char c;
      ifile.read(&c,1);
      int cval = (int)((unsigned char)c);
      pix->pixel_ref(x,y)=cval;
    }
  }
  return pix;
}


