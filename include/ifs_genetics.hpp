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

GenePoolRecordT genetical_optimize( size_t pool_size, 
				    size_t orphans_per_generation, 
				    size_t n_mutants, 
				    size_t n_crossovers,
				    PixelMap &sample, 
				    size_t generations,
				    size_t stop_if_no_improvement_after);

#endif
