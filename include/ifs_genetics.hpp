#ifndef __IFS_GENETICS_HPP_INCLUDED__
#define __IFS_GENETICS_HPP_INCLUDED__
#include <iostream>
#include "genetic_optim.hpp"

class Transform;
class PixelMap;
class Ruleset;

Transform merge_transforms(const Transform &t1, const Transform &t2, double p);
void normalize_pixmap( PixelMap &p );

class RulesetGenetics: public Genetics<Ruleset>
{
private:
  void mutate_global_noise(Ruleset &r);
  void mutate_insert(Ruleset &r);
  void mutate_delete(Ruleset &r);
  void mutate_modify(Ruleset &r);
  void random_modify_rule( Ruleset::Rule &r, double amount );

public:
  
  double noise_amount_global;
  double noise_amount_point;
  int crossover_size_jitter;
  size_t max_rules;

  RulesetGenetics();
  virtual Ruleset *orphan();
  virtual Ruleset *clone(const Ruleset &g);
  virtual Ruleset *mutant(const Ruleset &g);
  virtual Ruleset *crossover(const Ruleset &g1, const Ruleset &g2);
  void to_stream( std::ostream &os, const Ruleset&r );
  virtual void deallocate( Ruleset *g );
  
};

/**MEasures similarity of the image, produced by the IFS ruleset, and the given image,
   using cosine distance
   If image is normalized, then result if between 0 and 1.
*/
class CosineMeasureFitness: public FitnessFunction<Ruleset>{
  const PixelMap &sample;
  PixelMap canvas;
  point_t origin, size;
  double gamma;
public:
  size_t RENDER_STEPS_PER_PIXEL;
  CosineMeasureFitness( const PixelMap &sample_, const point_t &p0, const point_t &p1 );
  void set_gamma(double g){ gamma = g; };
  virtual double fitness(const Ruleset &rule);
};

#endif

