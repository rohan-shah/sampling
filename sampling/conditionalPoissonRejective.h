#ifndef CONDITIONAL_POISSON_REJECTIVE_HEADER_GUARD
#define CONDITIONAL_POISSON_REJECTIVE_HEADER_GUARD
#include <vector>
#include "includeMPFRSampling.h"
#include <boost/random/mersenne_twister.hpp>
#include "conditionalPoissonBase.h"
namespace sampling
{
	void conditionalPoissonRejective(conditionalPoissonArgs& args, boost::mt19937& randomSource);
}
#endif
