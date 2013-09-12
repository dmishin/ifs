#include "pgm.hpp"
#include <iostream>
#include <stdexcept>
#include <cmath>

size_t read_pgm_header_int(std::istream &ifile);
void skip_to_eol(std::istream &ifile);


void save_pgm( MonochromeImageReader &pixels, std::ostream &ofile )
{
  ofile << "P5\n"; //binary graymap
  ofile << pixels.width() << " " << pixels.height() << "\n";
  ofile << "255\n";

  const size_t BUF_SIZE = 1024;

  unsigned char buffer[BUF_SIZE];
  size_t sz = pixels.width()*pixels.height();
  size_t nWritten = 0;
  while( nWritten < sz ){
    size_t nRead = pixels.read_pixels(buffer, BUF_SIZE);
    ofile.write( (const char*)buffer, nRead );
    nWritten += nRead;
  };
}

void read_pgm_header( std::istream &ifile, size_t &w, size_t &h, size_t &maxcolor )
{
  char magic[2];
  ifile.read(&magic[0], 2);
  if (magic[0] !='P' || magic[1] != '5')
	  throw std::logic_error("Bad magic number: P5 was expected.");
  w = read_pgm_header_int(ifile);
  h = read_pgm_header_int(ifile);
  maxcolor = read_pgm_header_int(ifile);
}

void skip_to_eol(std::istream &ifile)
{
  char c;
  while(ifile){
    ifile.read( &c, 1);
    if (c=='\n' || c=='\r') return;
  }
}
size_t read_pgm_header_int(std::istream &ifile)
{
  char c;
  size_t x = 0;
  bool has_number = false;
  while(ifile){
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
void read_pgm( std::istream &ifile, MonochromeImageWriter &pixels )
{
  size_t w, h, maxcolor;
  read_pgm_header( ifile, w, h, maxcolor );
  if (maxcolor > 255)
    throw std::logic_error( "only 8-bit supported");
  if (w ==0 || h==0 )
    throw std::logic_error( "zero size of image");
  pixels.set_size(w,h);

  std::cerr<<"Reading image "<<w<<"x"<<h<<" with "<<maxcolor<<" colors\n";
  const size_t BUF_SIZE=1024;
  unsigned char buffer[BUF_SIZE];

  size_t sz = w*h;
  size_t nLeft = sz;
  while (ifile && (nLeft != 0)){
    ifile.read((char *)buffer, std::min(BUF_SIZE, nLeft));
    size_t nRead = ifile.gcount();
    pixels.write_pixels(buffer, nRead);
    nLeft -= nRead;
  }
  if (nLeft != 0)
	  throw std::logic_error("Premature end of file");
}


