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

const double GLOBAL_NOISE_AMOUNT = 0.05; //defines
const double POINT_NOISE_AMOUNT = 0.2; //defines
const int CROSSOVER_RANDOMIZE_SIZE = 2;
const size_t MAX_GENOME_SIZE = 10;
const size_t RENDER_STEPS_PER_PIXEL = 3;

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



void mutate_global_noise(Ruleset &r);
void mutate_insert(Ruleset &r);
void mutate_delete(Ruleset &r);
void mutate_modify(Ruleset &r);

void random_modify_rule( Ruleset::Rule &r, double amount )
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
void mutate_global_noise(Ruleset &r)
{
  using namespace std;
  double amount = GLOBAL_NOISE_AMOUNT * random_double();

  for (size_t i=0; i<r.rules.size(); ++i){
    random_modify_rule(r.rules[i], amount);
  }
  r.update_probabilities();
}
void mutate_insert(Ruleset &r)
{
  if (r.size() > MAX_GENOME_SIZE)
    mutate_delete(r);
  r.add( pow(random_double(), 3) );
  Transform &t(r.last().transform);
  t.offset = point_t(random_double()-0.5, random_double()-0.5);
  t.rot_scale( random_double()*3.1415926*2, random_double() );

  r.update_probabilities();
}

void mutate_delete(Ruleset &r)
{
  size_t max_size = 10;
  if (r.size() <= 2) return;
  if (r.size() >= max_size) return;
  size_t idx = rand() % r.size();
  r.rules.erase(r.rules.begin()+idx);
  r.update_probabilities();
}

void mutate_modify(Ruleset &r)
{
  double amount = POINT_NOISE_AMOUNT;
  if (r.size() ==0) return;
  size_t idx = rand() % r.size();
  random_modify_rule( r.rules[idx], amount );
  r.update_probabilities();
}

struct ByFitness{
  bool operator()(const GenePoolRecordT &r1, const GenePoolRecordT& r2)const{
    return r1.fitness > r2.fitness;
  }
};

std::ostream & operator << (std::ostream &os, const GenePoolRecordT &record)
{
  return
    os<<record.fitness
      <<" origin: "<<record.origin
      <<" born: "<<record.generation
      <<" sz:"<<record.genome->size();
}
GenePoolRecordT genetical_optimize( Genetics<Ruleset> &genetics,
				    size_t pool_size, 
				    size_t orphans_per_generation, 
				    size_t n_mutants, 
				    size_t n_crossovers,
				    FitnessFunction &fitness_func, 
				    size_t generations,
				    size_t stop_if_no_improvement_after)
{
  using namespace std;
  cout << "Starting genetical optimization with parameters:"<<endl;
  cout << "  pool size:"<<pool_size<<endl;
  cout << "  orphans per generation:"<<orphans_per_generation<<endl;
  cout << "  mutants per generation:"<<n_mutants<<endl;
  cout << "  crossovers per generation:"<<n_crossovers<<endl;

  GenePoolRecordT best;
  best.genome = NULL;
  best.fitness = -1;
    
  GenePoolT pool;
  for(size_t i=0; i<pool_size; ++i){
    pool.push_back( GenePoolRecordT( genetics.orphan(), "initial orphan") );
    pool.back().generation = 0;
  }

  for( size_t generation=0; generation < generations; ++generation){
    std::cout<<"Generation #"<<(generation+1)<<" of "<<generations<<std::endl;
    //add mutants    
    for(size_t i=0; i<n_mutants; ++i){
      size_t idx = rand()%pool_size;
      pool.push_back( GenePoolRecordT( genetics.mutant( *(pool[idx].genome) ),
				       "mutant") );
      pool.back().generation = generation;
    }
    //add crossovers
    for(size_t i=0; i<n_crossovers; ++i){
      size_t idx1 = rand()%pool_size;
      size_t idx2 = rand()%pool_size;
      pool.push_back( GenePoolRecordT( genetics.crossover( *(pool[idx1].genome),
							   *(pool[idx2].genome) ),
				       "crossover") );
      pool.back().generation = generation;
    }
    
    //add orphans
    for(size_t i=0; i<orphans_per_generation; ++i){
      pool.push_back( GenePoolRecordT( genetics.orphan(), "orphan" ) );
      pool.back().generation = generation;
    }
    //Update fitness for those who has not it.
    for( GenePoolT::iterator i=pool.begin(); i!=pool.end(); ++i){
      if (i->fitness >= 0) continue;//already calculated
      i->fitness = fitness_func.fitness( *(i->genome) );
    }
    //Update best sample
    if (pool.front().fitness > best.fitness || best.fitness < 0){
      best = pool.front();
      best.genome = genetics.clone(*best.genome);
    }
    //Remove the worst samples;
    std::sort(pool.begin(), pool.end(), ByFitness() );

    std::cout <<"   Best in history :"<<std::endl;
    std::cout<<"  "<<best<<std::endl;
    
    std::cout <<"   Best 3 fitness:"<<std::endl;
    for(size_t i=0; i<std::min((size_t)3,pool.size()); ++i){
      std::cout<<"  "<<i<<" "<<pool[i]<<std::endl;
    }

    if( pool.front().fitness < pool.back().fitness )
      throw std::logic_error("assertion: fitness sort is bad");

    //Filter old records
    
    if (true){
      size_t i=0;
      while(i < pool.size()){
	if (generation - pool[i].generation > 10){
	  genetics.deallocate(pool[i].genome);
	  pool.erase(pool.begin()+i);
	}else{
	  i += 1;
	}
      }
    }
    if (pool.size() > pool_size){
      for( GenePoolT::iterator i = pool.begin() + pool_size;
	   i < pool.end();
	   ++i ){
	genetics.deallocate(i->genome);
	i->genome = NULL;
	i->fitness = -1;
      }
      pool.erase(pool.begin() + pool_size, pool.end());
    }

    size_t youngest = 0;
    for(GenePoolT::iterator i=pool.begin(); i!= pool.end();++i){
      youngest = std::max(youngest, i->generation);
    }
    if (generation - youngest > stop_if_no_improvement_after){
      std::cout<<"Spent "<<stop_if_no_improvement_after<<" generations without improvement, stopping evolution"<<std::endl;
      break;
    }
  }
  for( GenePoolT::iterator i = pool.begin();
       i < pool.end();
       ++i ){
    delete (i->genome);
    i->genome = NULL;
    i->fitness = -1;
  }
  return best;
}

CosineMeasureFitness::CosineMeasureFitness( const PixelMap &sample_, const point_t &p0, const point_t &p1)
  :sample(sample_)
  ,canvas(sample.width, sample.height)
{
  origin = p0;
  size = p1 - p0;
  gamma = 5;
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
  size_t n = 2 + (rand()%(MAX_GENOME_SIZE - 1));
  for(size_t i=0; i<n; ++i){
    mutate_insert( *orphan );
  }
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
  return mutant;
}
Ruleset *RulesetGenetics::crossover(const Ruleset &r1, const Ruleset &r2)
{
  double merge_k = random_double()*2-0.5; //from -0.5 to 1.5
  Ruleset *crs = new Ruleset();
  
  int dn = rand()%(CROSSOVER_RANDOMIZE_SIZE*2+1) - CROSSOVER_RANDOMIZE_SIZE;
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
