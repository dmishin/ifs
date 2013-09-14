#ifndef __IFS_RAND_HPP_INCLUDED__
#define __IFS_RAND_HPP_INCLUDED__
#include "ruleset.hpp"

Transform merge_transforms(const Transform &t1, const Transform &t2, double p);

void normalize_pixmap( PixelMap &p );

#endif
