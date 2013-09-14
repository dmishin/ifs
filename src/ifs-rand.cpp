
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <time.h>

#include "ifs-rand.hpp"
#include "geometry.hpp"
#include "pgm.hpp"
#include "pixelmap.hpp"
#include "ruleset.hpp"
#include "ifs_genetics.hpp"

void render_sample_ruleset()
{
  srand( 11 ); //always the same
  Ruleset *r = make_orphan();

  PixelMap pix(800, 800);


  render_ruleset( pix, 
		  point_t(-1.5,-1.5),
		  point_t(3,3),
		  *r,
		  pix.width*pix.height*100 );

  pix.normalize(1);
  pix.apply_gamma(5);
  {
    std::ofstream out("sample-render-ruleset.pgm", std::ios::binary | std::ios::out);  
    PixelMapReader r(pix);
    save_pgm( r, out);
  }
}

int main( int argc, char *argv[] )
{
  srand((unsigned int)time(NULL));
  std::ifstream ifile("sample-small.pgm", std::ios::binary | std::ios::in);
  PixelMap pix1(0,0);
  {
    PixelMapWriter w(pix1);
    read_pgm(ifile, w);
  }
  normalize_pixmap( pix1 );
  CosineMeasureFitness fitness( pix1, point_t(-1,-1), point_t(1,1) );
  
  GenePoolRecordT result = 
    genetical_optimize( 10, //pool
			8, //orp
			15, //mut
			15, //cross
			fitness,
			30000,
			50);
  std::cout<<"Genetical optimization finished, rendering showing the result"<<std::endl;
  PixelMap pix2(800, 800);

  


  render_ruleset( pix2, 
		  point_t(-1.5,-1.5),
		  point_t(3,3),
		  *result.genome,
		  pix2.width*pix2.height*100 );
  pix2.normalize(1);
  pix2.apply_gamma(5);

  std::ofstream out("test-small.pgm", std::ios::binary | std::ios::out);  
  {
    PixelMapReader r(pix2);
    save_pgm( r, out);
  }

  //render_sample_ruleset();
  return 0;
}
