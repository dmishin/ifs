#include <cstdlib>
#include <algorithm>
#include <stdexcept>
#include "genetic_optim.hpp"
GenericGeneticalOptimizer::GenericGeneticalOptimizer(GenericGenetics &genetics_, GenericFitnessFunction &fitness_function_)
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
void GenericGeneticalOptimizer::set_parameters( size_t pool_size_,
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

void GenericGeneticalOptimizer::initialize_pool()
{
  for(size_t i=0; i<pool_size; ++i){
    pool.push_back( PoolRecord( genetics._orphan(), "initial orphan") );
    pool.back().generation = 0;
  }
}
void GenericGeneticalOptimizer::add_mutant()
{
  size_t idx = rand()%pool_size;
  pool.push_back( PoolRecord( genetics._mutant( pool[idx].genome ),
			      "mutant") );
  pool.back().generation = generation;
}
void GenericGeneticalOptimizer::add_crossover()
{
  size_t idx1 = rand()%pool_size;
  size_t idx2 = rand()%pool_size;
  pool.push_back( PoolRecord( genetics._crossover( pool[idx1].genome,
						  pool[idx2].genome ),
			      "crossover") );
  pool.back().generation = generation;
}
void GenericGeneticalOptimizer::add_orphan()
{
  pool.push_back( PoolRecord( genetics._orphan(), "orphan" ) );
  pool.back().generation = generation;
}
void GenericGeneticalOptimizer::update_fitness_values()
{
  for( PoolT::iterator i=pool.begin(); i!=pool.end(); ++i){
    if (i->fitness >= 0) continue;//already calculated
    i->fitness = fitness_function._fitness( i->genome );
  }
}

struct RecordsByFitness{
  bool operator()(const GenericGeneticalOptimizer::PoolRecord &r1, const GenericGeneticalOptimizer::PoolRecord& r2)const{
    return r1.fitness > r2.fitness;
  }
};
void GenericGeneticalOptimizer::run(size_t generations)
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
      best.genome = genetics._clone(best.genome);
    }

    //Remove the worst samples;

    std::cout <<"   Best in history :"<<std::endl
	      <<"  "<<show_record(best)<<std::endl
	      <<"   Best 3 fitness:"<<std::endl;

    for(size_t i=0; i<std::min((size_t)3,pool.size()); ++i){
      std::cout<<"  "<<i<<" "<<show_record(pool[i])<<std::endl;
    }

    if( pool.front().fitness < pool.back().fitness )
      throw std::logic_error("assertion: fitness sort is bad");

    //Filter old records
    
    if (die_if_older_than != 0){
      size_t i=0;
      while(i < pool.size()){
	if (generation - pool[i].generation > die_if_older_than){
	  genetics._deallocate(pool[i].genome);
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
	genetics._deallocate(i->genome);
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
void GenericGeneticalOptimizer::clear_pool()
{
  for( PoolT::iterator i = pool.begin();
       i < pool.end();
       ++i ){
    genetics._deallocate(i->genome);
    i->genome = NULL;
    i->fitness = -1;
  }
  pool.clear();
}

GenericGeneticalOptimizer::~GenericGeneticalOptimizer()
{
  clear_pool();
};



std::ostream & operator <<(std::ostream &os, const GenericGeneticalOptimizer::ShowRecord &show)
{
  os<<"fit: "<<show.record.fitness<<" origin:"<<show.record.origin<<" born:"<<show.record.generation<<" ";
  show.opt.genome_to_stream(os, show.record.genome);
  return os;
}
