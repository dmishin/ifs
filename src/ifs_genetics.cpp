#include <vector>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <stdlib.h>

#include "ifs-rand.hpp"
#include "geometry.hpp"
#include "pixelmap.hpp"
#include "ruleset.hpp"
#include "ifs_genetics.hpp"


Transform merge_transforms(const Transform &t1, const Transform &t2, double p)
{
  Transform t;
  for(size_t i=0; i<6; ++i){
    t.as_vector(i) = t1.as_vector(i)*p + t2.as_vector(i)*(1-p);
  }
  return t;
}


void normalize_pixmap( PixelMap &p )
{
  p.scale( 1.0 / p.norm2() );
}

double cosine_distance( const PixelMap &p1, const PixelMap &p2 )
{
  if (p1.width != p2.width) throw std::logic_error( "widths are different");
  if (p1.height != p2.height) throw std::logic_error( "heights are different");
  double s1=0, s12=0;
  size_t n = p1.pixels.size();
  for(size_t i=0; i<n; ++i){
    double v1 = p1.pixels[i];
    double v2 = p2.pixels[i];
    s1 += v1*v1;
    s12 += v1*v2;
  }
  if (s1 == 0.0) return 0;
  return s12 / sqrt(s1);
}




void RulesetGenetics::random_modify_rule( Ruleset::Rule &r, double amount )
{
    r.probability *= (1.0 + amount*(random_double()-0.5));
    Transform &t(r.transform); 

    t.offset += point_t(random_double()-0.5, random_double()-0.5) * amount;
    t.t00 *= (1.0 + amount*(random_double()-0.5));
    t.t01 *= (1.0 + amount*(random_double()-0.5));
    t.t10 *= (1.0 + amount*(random_double()-0.5));
    t.t11 *= (1.0 + amount*(random_double()-0.5));

    t.t00 += amount*(random_double()-0.5);
    t.t01 += amount*(random_double()-0.5);
    t.t10 += amount*(random_double()-0.5);
    t.t11 += amount*(random_double()-0.5);
}
void RulesetGenetics::mutate_global_noise(Ruleset &r)
{
  using namespace std;
  double amount = noise_amount_global * random_double();

  for (size_t i=0; i<r.rules.size(); ++i){
    random_modify_rule(r.rules[i], amount);
  }
}
void RulesetGenetics::mutate_insert(Ruleset &r)
{
  if (r.size() > max_rules)
    mutate_delete(r);
  r.add( pow(random_double(), 3) );
  Transform &t(r.last().transform);
  t.offset = point_t(random_double()-0.5, random_double()-0.5);
  t.rot_scale( random_double()*3.1415926*2, random_double() );

}

void RulesetGenetics::mutate_delete(Ruleset &r)
{
  size_t max_size = 10;
  if (r.size() <= 2) return;
  if (r.size() >= max_size) return;
  size_t idx = rand() % r.size();
  r.rules.erase(r.rules.begin()+idx);
}

void RulesetGenetics::mutate_modify(Ruleset &r)
{
  double amount = noise_amount_point;
  if (r.size() ==0) return;
  size_t idx = rand() % r.size();
  random_modify_rule( r.rules[idx], amount );
}

CosineMeasureFitness::CosineMeasureFitness( const PixelMap &sample_, const point_t &p0, const point_t &p1)
  :sample(sample_)
  ,canvas(sample.width, sample.height)
{
  origin = p0;
  size = p1 - p0;
  gamma = 5;
  RENDER_STEPS_PER_PIXEL = 3;
}

double CosineMeasureFitness::fitness(const Ruleset &rule)
{
  canvas.fill(0);
  render_ruleset( canvas,
		  origin, size,
		  rule,
		  canvas.width*canvas.height * RENDER_STEPS_PER_PIXEL );
  //pix.apply_treshold(RENDER_STEPS_PER_PIXEL * 0.5, 1);
  canvas.apply_gamma( gamma );
  return cosine_distance(canvas, sample);
}

Ruleset *RulesetGenetics::orphan()
{
  Ruleset * orphan = new Ruleset();
  size_t n = 2 + (rand()%(max_rules - 1));
  for(size_t i=0; i<n; ++i){
    mutate_insert( *orphan );
  }
  orphan->update_probabilities();
  return orphan;
}
Ruleset *RulesetGenetics::clone(const Ruleset &g)
{
  return new Ruleset(g);
}
Ruleset *RulesetGenetics::mutant(const Ruleset &g)
{
  Ruleset *mutant = clone(g);
  switch ( rand() % 4 ){
  case 0:
    mutate_global_noise( *mutant );
    break;
  case 1:
    mutate_insert(*mutant);
    break;
  case 2:
    mutate_delete(*mutant);
    break;
  case 3:
    mutate_modify(*mutant);
    break;
  default: 
    throw std::logic_error("Error: incorrect mutation type");
  }
  mutant->update_probabilities();

  return mutant;
}
Ruleset *RulesetGenetics::crossover(const Ruleset &r1, const Ruleset &r2)
{
  double merge_k = random_double()*2-0.5; //from -0.5 to 1.5
  Ruleset *crs = new Ruleset();
  
  int dn = rand()%(crossover_size_jitter*2+1) - crossover_size_jitter;
  int n1 = (r1.size() + r2.size()) / 2 + dn;
  size_t n = (n1 > 2)? n1 : 2;
  {
    std::back_insert_iterator<Ruleset::RulesT> inserter(crs->rules);
    std::copy(r1.rules.begin(), r1.rules.end(), inserter);
    std::copy(r2.rules.begin(), r2.rules.end(), inserter);
    for(Ruleset::RulesT::iterator i=crs->rules.begin(); i!=crs->rules.end();++i){
      i->probability *= 0.5;
    }
  }
  while( crs->size() > n ){
    size_t j, j1;
    crs->most_similar_rule_pair( j, j1 );

    if (j1 >= crs->size() ) break;
    double sp = crs->rules[j].probability + crs->rules[j1].probability;

    //weighted probabilities
    double wp1 = crs->rules[j].probability * merge_k;
    double wp2 = crs->rules[j1].probability * (1-merge_k);

    crs->rules[j].transform = merge_transforms(crs->rules[j].transform, 
					       crs->rules[j1].transform,
					       wp1 / (wp1 + wp2) );

    crs->rules[j].probability = sp;
    crs->rules.erase(crs->rules.begin() + j1 );
  }
  crs->update_probabilities();
  return crs;
}
void RulesetGenetics::deallocate( Ruleset *g )
{
  delete g;
};

RulesetGenetics::RulesetGenetics()
{
  noise_amount_global = 0.05; //defines
  noise_amount_point = 0.2; //defines
  crossover_size_jitter = 2;
  max_rules = 10;
}

struct RecordsByFitness{
  bool operator()(const GeneticalOptimizer::PoolRecord &r1, const GeneticalOptimizer::PoolRecord& r2)const{
    return r1.fitness > r2.fitness;
  }
};

GeneticalOptimizer::GeneticalOptimizer(Genetics<Ruleset> &genetics_, FitnessFunction<Ruleset> &fitness_function_)
  :genetics(genetics_)
  ,fitness_function(fitness_function_)
{
  set_parameters( 16,
		  5,
		  8,
		  8,
		  50,
		  10 );
  best.genome = NULL;
  best.fitness = -1;
  generation = 0;
}
void GeneticalOptimizer::set_parameters( size_t pool_size_,
					 size_t orphans_per_generation_, 
					 size_t n_mutants_, 
					 size_t n_crossovers_,
					 size_t stop_if_no_improvement_after_,
					 size_t die_if_older_than_)
{
  pool_size = pool_size_; 
  orphans_per_generation = orphans_per_generation_; 
  n_mutants = n_mutants_; 
  n_crossovers = n_crossovers_;
  stop_if_no_improvement_after = stop_if_no_improvement_after_;
  die_if_older_than = die_if_older_than_;
}

void GeneticalOptimizer::initialize_pool()
{
  for(size_t i=0; i<pool_size; ++i){
    pool.push_back( PoolRecord( genetics.orphan(), "initial orphan") );
    pool.back().generation = 0;
  }
}
void GeneticalOptimizer::add_mutant()
{
  size_t idx = rand()%pool_size;
  pool.push_back( PoolRecord( genetics.mutant( *(pool[idx].genome) ),
			      "mutant") );
  pool.back().generation = generation;
}
void GeneticalOptimizer::add_crossover()
{
  size_t idx1 = rand()%pool_size;
  size_t idx2 = rand()%pool_size;
  pool.push_back( PoolRecord( genetics.crossover( *(pool[idx1].genome),
						  *(pool[idx2].genome) ),
			      "crossover") );
  pool.back().generation = generation;
}
void GeneticalOptimizer::add_orphan()
{
  pool.push_back( PoolRecord( genetics.orphan(), "orphan" ) );
  pool.back().generation = generation;
}
void GeneticalOptimizer::update_fitness_values()
{
  for( PoolT::iterator i=pool.begin(); i!=pool.end(); ++i){
    if (i->fitness >= 0) continue;//already calculated
    i->fitness = fitness_function.fitness( *(i->genome) );
  }
}
std::ostream & operator << (std::ostream &os, const GeneticalOptimizer::PoolRecord &record)
{
  return
    os<<record.fitness
      <<" origin: "<<record.origin
      <<" born: "<<record.generation
      <<" sz:"<<record.genome->size();
}

void GeneticalOptimizer::run(size_t generations)
{
  using namespace std;
  cout << "  pool size:"<<pool_size<<endl
       << "  orphans per generation:"<<orphans_per_generation<<endl
       << "  mutants per generation:"<<n_mutants<<endl
       << "  crossovers per generation:"<<n_crossovers<<endl;

  if (pool.empty())
    initialize_pool();

  for( size_t iter=0; iter < generations; ++iter, ++generation){
    std::cout<<"Generation #"<<(generation+1)<<" of "<<generations<<std::endl;
    //add mutants    
    for(size_t i=0; i<n_mutants; ++i)
      add_mutant();

    //add crossovers
    for(size_t i=0; i<n_crossovers; ++i)
      add_crossover();
    
    //add orphans
    for(size_t i=0; i<orphans_per_generation; ++i)
      add_orphan();
    
    //Update fitness for those who has not it.
    update_fitness_values();

    sort(pool.begin(), pool.end(), RecordsByFitness() );

    //Update best sample
    if (pool.front().fitness > best.fitness || best.fitness < 0){
      best = pool.front();
      best.genome = genetics.clone(*best.genome);
    }

    //Remove the worst samples;

    std::cout <<"   Best in history :"<<std::endl
	      <<"  "<<best<<std::endl
	      <<"   Best 3 fitness:"<<std::endl;

    for(size_t i=0; i<std::min((size_t)3,pool.size()); ++i){
      std::cout<<"  "<<i<<" "<<pool[i]<<std::endl;
    }

    if( pool.front().fitness < pool.back().fitness )
      throw std::logic_error("assertion: fitness sort is bad");

    //Filter old records
    
    if (die_if_older_than != 0){
      size_t i=0;
      while(i < pool.size()){
	if (generation - pool[i].generation > die_if_older_than){
	  genetics.deallocate(pool[i].genome);
	  pool.erase(pool.begin()+i);
	}else{
	  i += 1;
	}
      }
    }
    if (pool.size() > pool_size){
      for( PoolT::iterator i = pool.begin() + pool_size;
	   i < pool.end();
	   ++i ){
	genetics.deallocate(i->genome);
	i->genome = NULL;
	i->fitness = -1;
      }
      pool.erase(pool.begin() + pool_size, pool.end());
    }

    size_t youngest = 0;
    for(PoolT::iterator i=pool.begin(); i!= pool.end();++i){
      youngest = max(youngest, i->generation);
    }
    if (generation - youngest > stop_if_no_improvement_after){
      cout<<"Spent "<<stop_if_no_improvement_after<<" generations without improvement, stopping evolution"<<std::endl;
      break;
    }
  }

}
void GeneticalOptimizer::clear_pool()
{
  for( PoolT::iterator i = pool.begin();
       i < pool.end();
       ++i ){
    delete (i->genome);
    i->genome = NULL;
    i->fitness = -1;
  }
  pool.clear();
}

GeneticalOptimizer::~GeneticalOptimizer()
{
  clear_pool();
};
