
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <cmath>
#include <cstdlib>

#include "geometry.hpp"
#include "pixelmap.hpp"
#include "ruleset.hpp"

bool Ruleset::most_similar_rule_pair( size_t &a, size_t &b)const
{
	if (rules.size() <= 1) return false;
	if (rules.size() == 2) {
		a = 0; b = 1;
		return true;
	}
	double dbest = -1;
	for(size_t i=0; i<rules.size(); ++i){
		for(size_t j=0; j<rules.size(); ++j){
			double d = rules[i].transform.distance(rules[j].transform);
			if (dbest < 0 || d < dbest){
				a = i;
				b = j;
				dbest = d;
			}
		}
	}
	return true;
}

size_t Ruleset::most_similar_rule( const Ruleset::Rule &r, size_t except_index )const
{
  size_t ibest = rules.size();
  double dbest = -1;
  for(size_t i=0; i<rules.size();++i){
    if (i==except_index) continue;
    double d = rules[i].transform.distance(r.transform);
    if (dbest < 0 || d < dbest){
      dbest = d;
      ibest = i;
    };
  };
  return ibest;
}
Ruleset::Rule & Ruleset::add( double dp )
{
  rules.push_back( Rule() );
  rules.back().probability = dp;
  return rules.back();
}
void Ruleset::update_probabilities()
{
  integral_probabilities.resize( rules.size() );
  if (rules.empty())
    return;

  double sp=0;
  for(size_t i=0; i<rules.size(); ++i){
    sp += rules[i].probability;
    integral_probabilities[i] = sp;
  }
  //now normalize the integral probabilities to make them sum to 1
  if (sp == 0.0) return;//avoid zero division.
  double inv_sp = 1.0 / sp;
  for(size_t i=0; i<integral_probabilities.size(); ++i){
    integral_probabilities[i] *= inv_sp;
  }
};

point_t Ruleset::apply( const point_t &p )const
{
  double r = random_double();
  //binary search in range [a;b)
  size_t a, b;
  a = 0; b = rules.size()-1;
  while (b - a > 1){ 
    //p(a) is always below r; p(b) is always above. 
    size_t c = (a+b)/2;
    if (integral_probabilities[c] <= r){
      a = c;
    }else{
      b = c-1;
    }
  };
  if (b-a == 1){
    if ( r <= integral_probabilities[a] )
      b=a;
  }
  return rules[b].transform.apply(p);
};

void Ruleset::apply_inplace( point_t &p )const
{
  double r = random_double();
  //binary search in range [a;b)
  size_t a, b;
  a = 0; b = rules.size()-1;
  while (b - a > 1){ 
    //p(a) is always below r; p(b) is always above. 
    size_t c = (a+b)/2;
    if (integral_probabilities[c] <= r){
      a = c;
    }else{
      b = c-1;
    }
  };
  if (b-a == 1){
    if ( r <= integral_probabilities[a] )
      b=a;
  }
  rules[b].transform.apply_inplace(p);
};


/*** //Bad idea. SLower than I expected
void render_ruleset_buf( PixelMap &pixels, 
		     const point_t &origin, 
		     const point_t &size,
		     const Ruleset &ruleset,
		     size_t n)
{
  const size_t BUF_SIZE = 1024*1024;
  point_t p(0,0);
  PixelMap::pixel_t * buffer[BUF_SIZE];
  size_t buffer_pos = 0;

  point_t scale( pixels.width / size.x, pixels.height / size.y );

  for(size_t i=0; i<n; ++i){
    p = ruleset.apply(p);
    if( fabs(p.x) > 1e3 || fabs(p.y) > 1e3 ){
      p = point_t(0,0);
    }
    point_t pp = (p-origin)*scale;
    int ix = (int)floor(pp.x);
    int iy = (int)floor(pp.y);
    if (pixels.contains(ix,iy)){
	  buffer[buffer_pos ++] = &pixels.pixel_ref(ix,iy);
	  if (buffer_pos >= BUF_SIZE){
		  PixelMap::pixel_t **pi=&buffer[0], **pend=buffer + buffer_pos;
		  std::sort(pi, pend);
		  for(;pi!=pend;++pi){
			  **pi += 1;
		  }
		  buffer_pos = 0;
	  }
	}
  }
  PixelMap::pixel_t **pi=buffer, **pend=buffer+buffer_pos;
  std::sort(pi, pend);
  for(;pi!=pend;++pi){
	**pi += 1;
  }
}
*/

void render_ruleset( PixelMap &pixels, 
			 const point_t &origin, 
			 const point_t &size,
			 const Ruleset &ruleset,
			 size_t n)
{
  point_t scale = point_t(pixels.width,pixels.height) / size;
  Transform logical_to_screen;
  logical_to_screen.offset = -origin * scale;
  logical_to_screen.set_scale( scale );
  point_t zero = logical_to_screen.apply(point_t(0,0)); //coordinates of the zero point in screen coordinate system

  //Build a ruleset, that works in the screen coordinates,
  //tu reduce number of coordinate transforms
  Ruleset tfm_ruleset;
  ruleset.transform_ruleset(logical_to_screen, tfm_ruleset);

  //start stochastic rendering
  point_t p = zero;
  for(size_t i=0; i<n; ++i){
    tfm_ruleset.apply_inplace(p);

    if( fabs(p.x) > 1e6 || fabs(p.y) > 1e6 ){
      p = zero; 
    }
    int ix = (int)(p.x); //without floor it is faster
    int iy = (int)(p.y);
    if (pixels.contains(ix,iy)){
      pixels.pixel_ref(ix,iy) += 1;
    }
  }
  
}

//T is transform from the old coordinate system to the new.
void Ruleset::transform_ruleset( const Transform &t, Ruleset &out )const
{
  Transform inv_t = t.inverse();
  out.integral_probabilities = integral_probabilities;
  out.rules = rules;
  for(RulesT::iterator i=out.rules.begin(), e=out.rules.end(); i!=e; ++i ){
    i->transform = t.apply(i->transform.apply(inv_t));
  }
}
