#ifndef __IFS_GENETICS_HPP_INCLUDED__
#define __IFS_GENETICS_HPP_INCLUDED__
#include <iostream>

typedef void *GenericGenomePtr;
typedef const void *GenericGenomeCPtr;

class GenericFitnessFunction{
protected:
	virtual double _fitness(GenericGenomeCPtr g)=0;
	template <typename T> friend class GeneticalOptimizer;
	friend class GenericGeneticalOptimizer;
};

/**FItness function evaluates genome and determines, how good it is.
   It must return positive value;
*/   
template< typename Genome >
class FitnessFunction: public GenericFitnessFunction{
protected:
  virtual double _fitness(GenericGenomeCPtr g){return fitness(*(Genome*)g);};
public:
  virtual double fitness(const Genome &rule)=0;
  friend class GenericGeneticalOptimizer;
};
class GenericGenetics{
protected:
  virtual GenericGenomePtr _orphan()=0;
  virtual GenericGenomePtr _clone(GenericGenomeCPtr g)=0;
  virtual GenericGenomePtr _mutant(GenericGenomeCPtr g)=0;
  virtual GenericGenomePtr _crossover(GenericGenomeCPtr g1, GenericGenomeCPtr g2)=0;
  virtual void _deallocate( GenericGenomePtr )=0;
  friend class GenericGeneticalOptimizer;
  friend class ShowGenericGenome;
};

/**Genetics is collection of operations over genomes.
 */
template< typename Genome >
class Genetics: public GenericGenetics{
private:
  const Genome *from_generic(GenericGenomeCPtr p)const{ return (const Genome*)p; };
  Genome *from_generic(GenericGenomePtr p)const{ return (Genome*)p; };
  GenericGenomePtr to_generic( Genome *p )const{ return (GenericGenomePtr)p; };

protected:

  virtual GenericGenomePtr _orphan()
    {return to_generic(orphan()); };
  virtual GenericGenomePtr _clone(GenericGenomeCPtr g)
   {return to_generic(clone( *from_generic(g)));};
  virtual GenericGenomePtr _mutant(GenericGenomeCPtr g)
   {return to_generic(mutant( *from_generic(g) ));};
  virtual GenericGenomePtr _crossover(GenericGenomeCPtr g1, GenericGenomeCPtr g2)
  {return to_generic(crossover(*from_generic(g1), *from_generic(g2)));};
  virtual void _deallocate( GenericGenomePtr g)
    {deallocate( from_generic(g));};
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

class GenericGeneticalOptimizer{
public:
  struct PoolRecord{
    GenericGenomePtr genome;
    double fitness;
    std::string origin;
    size_t generation;
    PoolRecord():genome(NULL), fitness(-1){};
    PoolRecord(GenericGenomePtr g): genome(g), fitness(-1){};
    PoolRecord(GenericGenomePtr g, const std::string &o): genome(g), fitness(-1), origin(o){};
  };

  struct ShowRecord{
    const GenericGeneticalOptimizer &opt;
    const PoolRecord &record;
    ShowRecord(const GenericGeneticalOptimizer &o, const PoolRecord &r):opt(o),record(r){};
  };
  ShowRecord show_record(const PoolRecord &r)const{ return ShowRecord(*this, r); };

  typedef std::vector<PoolRecord> PoolT;
private:
  size_t pool_size; 
  size_t orphans_per_generation; 
  size_t n_mutants; 
  size_t n_crossovers;
  size_t generations;
  size_t stop_if_no_improvement_after;
  size_t die_if_older_than;

  GenericGenetics &genetics;
  GenericFitnessFunction &fitness_function;
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
protected:
  virtual void genome_to_stream(std::ostream &os, GenericGenomeCPtr g)const=0;
public:
  GenericGeneticalOptimizer(GenericGenetics &genetics_, GenericFitnessFunction &fitness_function_);
  ~GenericGeneticalOptimizer();
  void set_parameters( size_t pool_size_,
		       size_t orphans_per_generation_, 
		       size_t n_mutants_, 
		       size_t n_crossovers_,
		       size_t stop_if_no_improvement_after_,
		       size_t die_if_older_than_);
  void run(size_t generations);
  const PoolRecord &get_best()const{ return best; };

  friend std::ostream & operator <<(std::ostream &os, const GenericGeneticalOptimizer::ShowRecord &show);
};

std::ostream & operator <<(std::ostream &os, const GenericGeneticalOptimizer::ShowRecord &show);

template< typename Genome >
class GeneticalOptimizer: public GenericGeneticalOptimizer{
protected:
  virtual void genome_to_stream(std::ostream &os, GenericGenomeCPtr g)const{ genome_to_stream(os, *(Genome*)g); };
public:
  virtual void genome_to_stream(std::ostream &os, const Genome &g)const{ os<<"[genome]"; };
  GeneticalOptimizer(Genetics<Genome> &g, FitnessFunction<Genome> &f)
    :GenericGeneticalOptimizer(g,f){};
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

