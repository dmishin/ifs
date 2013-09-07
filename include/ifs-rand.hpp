#ifndef __IFS_RAND_HPP_INCLUDED__
#define __IFS_RAND_HPP_INCLUDED__
#include <vector>
#include "geometry.hpp"
#include "util.hpp"

Transform merge_transforms(const Transform &t1, const Transform &t2, double p);

/**Describes fractal*/
class Ruleset{
public:
  class Rule{
  public:
    double probability;
    Transform transform;
  };
  std::vector<double> integral_probabilities;
  typedef std::vector<Rule> RulesT;
  RulesT rules;
  void update_probabilities();
  double last_p()const;
  Rule & add( double dp );
  point_t apply( const point_t p )const;
  Rule &last(){ return rules.back(); };
  const Rule &last()const{ return rules.back(); };
  size_t size()const{ return rules.size(); };
  size_t most_similar_rule( const Ruleset::Rule &r, size_t except_index )const;
};


class PixelMap{
public:
  typedef float pixel_t;
  std::vector<pixel_t> pixels;
  size_t width, height;

  PixelMap( size_t w, size_t h );
  void fill( pixel_t v );

  size_t pixel_idx(size_t x, size_t y)const{ return x+y*width; };
  pixel_t &pixel_ref(int x, int y){ return pixels[pixel_idx(x,y)]; };
  bool contains( int x, int y )const;
  pixel_t max_value()const;
  void scale( double k );
};

void render_ruleset( PixelMap &pixels, 
		     const point_t &origin, 
		     const point_t &size,
		     const Ruleset &ruleset,
		     size_t n);

void normalize_pixmap( PixelMap &p );

#endif
