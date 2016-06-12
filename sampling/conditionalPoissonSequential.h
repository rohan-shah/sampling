#ifndef CONDITIONAL_POISSON_SEQUENTIAL_HEADER_GUARD
#define CONDITIONAL_POISSON_SEQUENTIAL_HEADER_GUARD
#include "conditionalPoissonBase.h"
#include <boost/random/mersenne_twister.hpp>
namespace sampling
{
	struct conditionalPoissonSequentialArgs : public conditionalPoissonArgs
	{
		std::vector<mpfr_class>* inclusionProbabilities;
	};
	void conditionalPoissonSequential(conditionalPoissonSequentialArgs& args, boost::mt19937& randomSource);
}
#endif
