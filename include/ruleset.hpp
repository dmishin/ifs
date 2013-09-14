#ifndef __RULESET_HPP_INCLDUED__
#define __RULESET_HPP_INCLDUED__

#include <vector>
#include "geometry.hpp"
#include "util.hpp"

/**Describes fractal*/
class Ruleset{
public:
  struct Rule{
    double probability;
    Transform transform;
  };

  std::vector<double> integral_probabilities;
  typedef std::vector<Rule> RulesT;
  RulesT rules;

public:
  void update_probabilities();
  Rule & add( double dp );
  point_t apply( const point_t &p )const;
  void apply_inplace( point_t &p )const;

  Rule &last(){ return rules.back(); };
  const Rule &last()const{ return rules.back(); };
  size_t size()const{ return rules.size(); };
  size_t most_similar_rule( const Ruleset::Rule &r, size_t except_index )const;
  bool most_similar_rule_pair( size_t &a, size_t &b)const;
  void transform_ruleset( const Transform &t, Ruleset &out )const;
};


class PixelMap;
void render_ruleset( PixelMap &pixels, 
		     const point_t &origin, 
		     const point_t &size,
		     const Ruleset &ruleset,
		     size_t n);


#endif
