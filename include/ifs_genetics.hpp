#ifndef __IFS_GENETICS_HPP_INCLUDED__
#define __IFS_GENETICS_HPP_INCLUDED__

class Transform;
class PixelMap;
class Ruleset;

Transform merge_transforms(const Transform &t1, const Transform &t2, double p);
void normalize_pixmap( PixelMap &p );

Ruleset * make_clone( const Ruleset &r );
Ruleset * make_orphan();
Ruleset * mutate( const Ruleset &r );
Ruleset * crossover( const Ruleset &r1, const Ruleset &r2);


struct GenePoolRecordT{
  Ruleset *genome;
  double fitness;
  std::string origin;
  size_t generation;
  GenePoolRecordT():genome(NULL), fitness(-1){};
  GenePoolRecordT(Ruleset *g): genome(g), fitness(-1){};
  GenePoolRecordT(Ruleset *g, const std::string &o): genome(g), fitness(-1), origin(o){};
};

typedef std::vector<GenePoolRecordT> GenePoolT;
class FitnessFunction;
GenePoolRecordT genetical_optimize( size_t pool_size, 
				    size_t orphans_per_generation, 
				    size_t n_mutants, 
				    size_t n_crossovers,
				    FitnessFunction &fitness_function, 
				    size_t generations,
				    size_t stop_if_no_improvement_after);

class FitnessFunction{
public:
  virtual double fitness(const Ruleset &rule)=0;
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
