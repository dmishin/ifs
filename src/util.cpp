#include <stdlib.h>
#include "util.hpp"
double random_double()
{
  return rand() / (double) RAND_MAX;
}
