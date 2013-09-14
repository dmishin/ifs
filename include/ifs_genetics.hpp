#ifndef __IFS_GENETICS_HPP_INCLUDED__
#define __IFS_GENETICS_HPP_INCLUDED__
#include <iostream>

/**FItness function evaluates genome and determines, how good it is.
   It must return positive value;
*/   
template< typename Genome >
class FitnessFunction{
public:
  virtual double fitness(const Genome &rule)=0;
};

/**Genetics is collection of operations over genomes.
 */
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

GenePoolRecordT genetical_optimize( Genetics<Ruleset> &genetics,
				    size_t pool_size, 
				    size_t orphans_per_generation, 
				    size_t n_mutants, 
				    size_t n_crossovers,
				    FitnessFunction<Ruleset> &fitness_function, 
				    size_t generations,
				    size_t stop_if_no_improvement_after);

class GeneticalOptimizer{
public:
  struct PoolRecord{
    Ruleset *genome;
    double fitness;
    std::string origin;
    size_t generation;
    PoolRecord():genome(NULL), fitness(-1){};
    PoolRecord(Ruleset *g): genome(g), fitness(-1){};
    PoolRecord(Ruleset *g, const std::string &o): genome(g), fitness(-1), origin(o){};
  };
  typedef std::vector<PoolRecord> PoolT;
private:
  size_t pool_size; 
  size_t orphans_per_generation; 
  size_t n_mutants; 
  size_t n_crossovers;
  size_t generations;
  size_t stop_if_no_improvement_after;
  size_t die_if_older_than;

  Genetics<Ruleset> &genetics;
  FitnessFunction<Ruleset> &fitness_function;
  PoolT pool;
  PoolRecord best;
  size_t generation;
private:
  void initialize_pool();
  void add_mutant();
  void add_crossover();
  void add_orphan();
  void update_fitness_values();
  void clear_pool();
public:
  GeneticalOptimizer(Genetics<Ruleset> &genetics_, FitnessFunction<Ruleset> &fitness_function_);
  void set_parameters( size_t pool_size_,
		       size_t orphans_per_generation_, 
		       size_t n_mutants_, 
		       size_t n_crossovers_,
		       size_t stop_if_no_improvement_after_,
		       size_t die_if_older_than_);
  void run(size_t generations);
};
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
