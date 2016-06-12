#ifndef CONDITIONAL_POISSON_DRAFTING_HEADER_GUARD
#define CONDITIONAL_POISSON_DRAFTING_HEADER_GUARD
#include <vector>
#include "includeMPFRSampling.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include "conditionalPoissonBase.h"
namespace sampling
{
	struct conditionalPoissonDraftingArgs : public conditionalPoissonArgs
	{
		std::vector<mpfr_class>* inclusionProbabilities;
		//This extra member is needed for recomputing the drafting probabilities after a unit is selected. 
		std::vector<mpfr_class> inclusionProbabilities2, inclusionProbabilities3;
	};
	void conditionalPoissonDrafting(conditionalPoissonDraftingArgs& args, boost::mt19937& randomSource);
}
#endif
