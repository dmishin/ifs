#ifndef __IFS_GENETICS_HPP_INCLUDED__
#define __IFS_GENETICS_HPP_INCLUDED__
#include <iostream>

class FitnessFunction{
public:
  virtual double fitness(const Ruleset &rule)=0;
};
template< typename Genome >
class Genetics{
public:
  virtual Genome *orphan()=0;
  virtual Genome *clone(const Genome &g)=0;
  virtual Genome *mutant(const Genome &g)=0;
  virtual Genome *crossover(const Genome &g1, const Genome &g2)=0;
  virtual void deallocate( Genome *g )=0;
};

class Transform;
class PixelMap;
class Ruleset;

Transform merge_transforms(const Transform &t1, const Transform &t2, double p);
void normalize_pixmap( PixelMap &p );

struct GenePoolRecordT{
  Ruleset *genome;
  double fitness;
  std::string origin;
  size_t generation;
  GenePoolRecordT():genome(NULL), fitness(-1){};
  GenePoolRecordT(Ruleset *g): genome(g), fitness(-1){};
  GenePoolRecordT(Ruleset *g, const std::string &o): genome(g), fitness(-1), origin(o){};
};

std::ostream & operator << (std::ostream &os, const GenePoolRecordT &record);

typedef std::vector<GenePoolRecordT> GenePoolT;
class FitnessFunction;
GenePoolRecordT genetical_optimize( Genetics<Ruleset> &genetics,
				    size_t pool_size, 
				    size_t orphans_per_generation, 
				    size_t n_mutants, 
				    size_t n_crossovers,
				    FitnessFunction &fitness_function, 
				    size_t generations,
				    size_t stop_if_no_improvement_after);


class RulesetGenetics: public Genetics<Ruleset>
{
public:
  virtual Ruleset *orphan();
  virtual Ruleset *clone(const Ruleset &g);
  virtual Ruleset *mutant(const Ruleset &g);
  virtual Ruleset *crossover(const Ruleset &g1, const Ruleset &g2);
  virtual void deallocate( Ruleset *g );
  
};
class CosineMeasureFitness: public FitnessFunction{
  const PixelMap &sample;
  PixelMap canvas;
  point_t origin, size;
  double gamma;
public:
  CosineMeasureFitness( const PixelMap &sample_, const point_t &p0, const point_t &p1 );
  void set_gamma(double g){ gamma = g; };
  virtual double fitness(const Ruleset &rule);
};

#endif
