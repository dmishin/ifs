
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <time.h>

#include "ifs-rand.hpp"
#include "geometry.hpp"
#include "pgm.hpp"
#include "pixelmap.hpp"

const double GLOBAL_NOISE_AMOUNT = 0.005; //defines
const double POINT_NOISE_AMOUNT = 0.1; //defines
const int CROSSOVER_RANDOMIZE_SIZE = 2;
const size_t MAX_GENOME_SIZE = 7;
const size_t RENDER_STEPS_PER_PIXEL = 2;

bool Ruleset::most_similar_rule_pair( size_t &a, size_t &b)const
{
	if (rules.size() <= 1) return false;
	if (rules.size() == 2) {
		a = 0; b = 1;
		return true;
	}
	double dbest = -1;
	for(size_t i=0; i<rules.size(); ++i){
		for(size_t j=0; j<rules.size(); ++j){
			double d = rules[i].transform.distance(rules[j].transform);
			if (dbest < 0 || d < dbest){
				a = i;
				b = j;
				dbest = d;
			}
		}
	}
	return true;
}

size_t Ruleset::most_similar_rule( const Ruleset::Rule &r, size_t except_index )const
{
  size_t ibest = rules.size();
  double dbest = -1;
  for(size_t i=0; i<rules.size();++i){
    if (i==except_index) continue;
    double d = rules[i].transform.distance(r.transform);
    if (dbest < 0 || d < dbest){
      dbest = d;
      ibest = i;
    };
  };
  return ibest;
}
Transform merge_transforms(const Transform &t1, const Transform &t2, double p)
{
  Transform t;
  for(size_t i=0; i<6; ++i){
    t.as_vector(i) = t1.as_vector(i)*p + t2.as_vector(i)*(1-p);
  }
  return t;
}
double Ruleset::last_p()const
{
  if (integral_probabilities.empty()) return 0;
  return integral_probabilities.back();
}

Ruleset::Rule & Ruleset::add( double dp )
{
  double ip = last_p()+dp;
  rules.push_back( Rule() );
  rules.back().probability = dp;
  integral_probabilities.push_back(ip);
  return rules.back();
}
void Ruleset::update_probabilities()
{
  integral_probabilities.resize( rules.size() );
  double sp=0;
  for(size_t i=0; i<rules.size(); ++i){
    sp += rules[i].probability;
    integral_probabilities[i] = sp;
  }
};

point_t Ruleset::apply( const point_t &p )const
{
  double r = random_double() * last_p();
  //binary search in range [a;b)
  size_t a, b;
  a = 0; b = rules.size()-1;
  while (b - a > 1){ 
    //p(a) is always below r; p(b) is always above. 
    size_t c = (a+b)/2;
    if (integral_probabilities[c] <= r){
      a = c;
    }else{
      b = c-1;
    }
  };
  if (b-a == 1){
    if ( r <= integral_probabilities[a] )
      b=a;
  }
  return rules[b].transform.apply(p);
};

void Ruleset::apply_inplace( point_t &p )const
{
  double r = random_double() * last_p();
  //binary search in range [a;b)
  size_t a, b;
  a = 0; b = rules.size()-1;
  while (b - a > 1){ 
    //p(a) is always below r; p(b) is always above. 
    size_t c = (a+b)/2;
    if (integral_probabilities[c] <= r){
      a = c;
    }else{
      b = c-1;
    }
  };
  if (b-a == 1){
    if ( r <= integral_probabilities[a] )
      b=a;
  }
  rules[b].transform.apply_inplace(p);
};


/*** //Bad idea. SLower than I expected
void render_ruleset_buf( PixelMap &pixels, 
		     const point_t &origin, 
		     const point_t &size,
		     const Ruleset &ruleset,
		     size_t n)
{
  const size_t BUF_SIZE = 1024*1024;
  point_t p(0,0);
  PixelMap::pixel_t * buffer[BUF_SIZE];
  size_t buffer_pos = 0;

  point_t scale( pixels.width / size.x, pixels.height / size.y );

  for(size_t i=0; i<n; ++i){
    p = ruleset.apply(p);
    if( fabs(p.x) > 1e3 || fabs(p.y) > 1e3 ){
      p = point_t(0,0);
    }
    point_t pp = (p-origin)*scale;
    int ix = (int)floor(pp.x);
    int iy = (int)floor(pp.y);
    if (pixels.contains(ix,iy)){
	  buffer[buffer_pos ++] = &pixels.pixel_ref(ix,iy);
	  if (buffer_pos >= BUF_SIZE){
		  PixelMap::pixel_t **pi=&buffer[0], **pend=buffer + buffer_pos;
		  std::sort(pi, pend);
		  for(;pi!=pend;++pi){
			  **pi += 1;
		  }
		  buffer_pos = 0;
	  }
	}
  }
  PixelMap::pixel_t **pi=buffer, **pend=buffer+buffer_pos;
  std::sort(pi, pend);
  for(;pi!=pend;++pi){
	**pi += 1;
  }
}
*/

void render_ruleset( PixelMap &pixels, 
			 const point_t &origin, 
			 const point_t &size,
			 const Ruleset &ruleset,
			 size_t n)
{
  point_t scale = point_t(pixels.width,pixels.height) / size;
  Transform logical_to_screen;
  logical_to_screen.offset = -origin * scale;
  logical_to_screen.set_scale( scale );
  point_t zero = logical_to_screen.apply(point_t(0,0)); //coordinates of the zero point in screen coordinate system

  //Build a ruleset, that works in the screen coordinates,
  //tu reduce number of coordinate transforms
  Ruleset tfm_ruleset;
  ruleset.transform_ruleset(logical_to_screen, tfm_ruleset);

  //start stochastic rendering
  point_t p = zero;
  for(size_t i=0; i<n; ++i){
    tfm_ruleset.apply_inplace(p);

    if( fabs(p.x) > 1e6 || fabs(p.y) > 1e6 ){
      p = zero; 
    }
    int ix = (int)(p.x); //without floor it is faster
    int iy = (int)(p.y);
    if (pixels.contains(ix,iy)){
      pixels.pixel_ref(ix,iy) += 1;
    }
  }
  
}

//T is transform from the old coordinate system to the new.
void Ruleset::transform_ruleset( const Transform &t, Ruleset &out )const
{
  Transform inv_t = t.inverse();
  out.integral_probabilities = integral_probabilities;
  out.rules = rules;
  for(RulesT::iterator i=out.rules.begin(), e=out.rules.end(); i!=e; ++i ){
    i->transform = t.apply(i->transform.apply(inv_t));
  }
}


void normalize_pixmap( PixelMap &p )
{
  double s=0;
  size_t n = p.pixels.size();
  for(size_t i=0; i<n; ++i){
    s += sqr(p.pixels[i]);
  }
  p.scale( 1.0 / sqrt(s) );
}

double angle_measure( const PixelMap &p1, const PixelMap &p2 )
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

Ruleset * make_orphan()
{
  Ruleset * orphan = new Ruleset();
  size_t n = 2 + (rand()%(MAX_GENOME_SIZE - 1));
  for(size_t i=0; i<n; ++i){
    mutate_insert( *orphan );
  }
  return orphan;
}

Ruleset * mutate( const Ruleset &r )
{
  Ruleset * clone = new Ruleset(r);

  switch ( rand() % 4 ){
  case 0:
    mutate_global_noise( *clone );
    break;
  case 1:
    mutate_insert(*clone);
    break;
  case 2:
    mutate_delete(*clone);
    break;
  case 3:
    mutate_modify(*clone);
    break;
  default: 
    throw std::logic_error("Error: incorrect mutation type");
  }
  return clone;
}
template<typename T>
T cap( T a, T b, T x )
{
  using namespace std;
  return min(max(x, a), b);
}

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
}
void mutate_modify(Ruleset &r)
{
  double amount = POINT_NOISE_AMOUNT;
  if (r.size() ==0) return;
  size_t idx = rand() % r.size();
  random_modify_rule( r.rules[idx], amount );
  r.update_probabilities();
}

Ruleset * crossover( const Ruleset &r1, const Ruleset &r2)
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

	double wp1 = crs->rules[j].probability * merge_k;
	double wp2 = crs->rules[j1].probability * (1-merge_k);

    crs->rules[j].transform = merge_transforms( 
		crs->rules[j].transform, 
		crs->rules[j1].transform,
		wp1 / (wp1 + wp2) );

    crs->rules[j].probability = sp;
    crs->rules.erase(crs->rules.begin() + j1 );
  }
  crs->update_probabilities();
  return crs;
}
struct GenePoolRecordT{
  Ruleset *genome;
  double fitness;
  std::string origin;
  size_t generation;
  GenePoolRecordT():genome(NULL), fitness(-1){};
  GenePoolRecordT(Ruleset *g): genome(g), fitness(-1){};
  GenePoolRecordT(Ruleset *g, const std::string &o): genome(g), fitness(-1), origin(o){};
};
struct ByFitness{
  bool operator()(const GenePoolRecordT &r1, const GenePoolRecordT& r2)const{
    return r1.fitness > r2.fitness;
  }
};

typedef std::vector<GenePoolRecordT> GenePoolT;

GenePoolRecordT genetical_optimize( size_t pool_size, 
				    size_t orphans_per_generation, 
				    size_t n_mutants, 
				    size_t n_crossovers,
				    PixelMap &sample, 
				    size_t generations,
				    size_t stop_if_no_improvement_after)
{
  using namespace std;
  cout << "Starting genetical optimization with parameters:"<<endl;
  cout << "  pool size:"<<pool_size<<endl;
  cout << "  orphans per generation:"<<orphans_per_generation<<endl;
  cout << "  mutants per generation:"<<n_mutants<<endl;
  cout << "  crossovers per generation:"<<n_crossovers<<endl;
  PixelMap pix(sample.width, sample.height);

  GenePoolT pool;
  for(size_t i=0; i<pool_size; ++i){
    pool.push_back( GenePoolRecordT( make_orphan(), "initial orphan") );
    pool.back().generation = 0;
  }

  for( size_t generation=0; generation < generations; ++generation){
    std::cout<<"Generation #"<<(generation+1)<<" of "<<generations<<std::endl;
    //add mutants    
    for(size_t i=0; i<n_mutants; ++i){
      size_t idx = rand()%pool_size;
      pool.push_back( GenePoolRecordT( mutate( *(pool[idx].genome) ),
				       "mutant") );
      pool.back().generation = generation;
    }
    //add crossovers
    for(size_t i=0; i<n_crossovers; ++i){
      size_t idx1 = rand()%pool_size;
      size_t idx2 = rand()%pool_size;
      pool.push_back( GenePoolRecordT( crossover( *(pool[idx1].genome),
						  *(pool[idx2].genome) ),
				       "crossover") );
      pool.back().generation = generation;
    }
    
    //add orphans
    for(size_t i=0; i<orphans_per_generation; ++i){
      pool.push_back( GenePoolRecordT( make_orphan(), "orphan" ) );
      pool.back().generation = generation;
    }
    //Update fitness for those who has not it.
    for( GenePoolT::iterator i=pool.begin(); i!=pool.end(); ++i){
      if (i->fitness >= 0) continue;//already calculated
      pix.fill(0);
      render_ruleset( pix, 
		      point_t(-1,-1),
		      point_t(2,2),
		      *(i->genome),
		      pix.width*pix.height*RENDER_STEPS_PER_PIXEL );
	  //pix.apply_treshold(RENDER_STEPS_PER_PIXEL * 0.5, 1);
	  pix.apply_gamma(5);
      i->fitness = angle_measure(pix, sample);
    }
    //Remove the worst samples;
    std::sort(pool.begin(), pool.end(), ByFitness() );

    std::cout <<"   Best 3 fitness:"<<std::endl;
    for(size_t i=0; i<std::min((size_t)3,pool.size()); ++i)
      std::cout<<"       "
	       <<i<<" "<<pool[i].fitness
	       <<" origin: "<<pool[i].origin
	       <<" born at: "<<pool[i].generation
	       <<" sz:"<<pool[i].genome->size()
	       <<std::endl;

    if( pool.front().fitness < pool.back().fitness )
      throw std::logic_error("assertion: fitness sort is bad");

    //Filter old records
    
    if (true){
      size_t i=0;
      while(i < pool.size()){
		  if (generation - pool[i].generation > 10){
			  delete pool[i].genome;
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
	delete (i->genome);
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
  for( GenePoolT::iterator i = pool.begin() + 1;
       i < pool.end();
       ++i ){
    delete (i->genome);
    i->genome = NULL;
    i->fitness = -1;
  }
  return pool[0];
}

void render_sample_ruleset()
{
  srand( 11 ); //always the same
  Ruleset *r = make_orphan();

  PixelMap pix(800, 800);


  render_ruleset( pix, 
		  point_t(-1.5,-1.5),
		  point_t(3,3),
		  *r,
		  pix.width*pix.height*100 );

  pix.normalize(1);
  pix.apply_gamma(5);
  {
    std::ofstream out("sample-render-ruleset.pgm", std::ios::binary | std::ios::out);  
    PixelMapReader r(pix);
    save_pgm( r, out);
  }
}
int main( int argc, char *argv[] )
{
  /*
  srand((unsigned int)time(NULL));
  std::ifstream ifile("sample-small.pgm", std::ios::binary | std::ios::in);
  PixelMap pix1(0,0);
  {
    PixelMapWriter w(pix1);
    read_pgm(ifile, w);
  }
  normalize_pixmap( pix1 );
  
  GenePoolRecordT result = 
    genetical_optimize( 10, //pool
			8, //orp
			1, //mut
			32, //cross
			pix1,
			300,
			50);
  std::cout<<"Genetical optimization finished, rendering showing the result"<<std::endl;
  PixelMap pix2(800, 800);


  render_ruleset( pix2, 
		  point_t(-1.5,-1.5),
		  point_t(3,3),
		  ruleset,
		  pix2.width*pix2.height*100 );
  pix2.normalize(1);
  pix2.apply_gamma(5);

  std::ofstream out("test-small.pgm", std::ios::binary | std::ios::out);  
  {
    PixelMapReader r(pix2);
    save_pgm( r, out);
  }
  */

  render_sample_ruleset();
  return 0;
}
