#ifndef __PGM__HPP_INCLUDED__
#define __PGM__HPP_INCLUDED__
#include <iostream>

class PixelMap;

void save_pgm( const PixelMap &pixels, std::ostream &ofile, double gamma=1.0 );
PixelMap * read_pgm( std::istream &ifile );

#endif
