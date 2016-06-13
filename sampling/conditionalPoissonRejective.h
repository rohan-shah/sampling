#ifndef CONDITIONAL_POISSON_REJECTIVE_HEADER_GUARD
#define CONDITIONAL_POISSON_REJECTIVE_HEADER_GUARD
#include <vector>
#include "includeMPFRSampling.h"
#include <boost/random/mersenne_twister.hpp>
#include "conditionalPoissonBase.h"
namespace sampling
{
	struct conditionalPoissonRejectiveArgs : public conditionalPoissonArgs
	{
	public:
		conditionalPoissonRejectiveArgs(bool calculateInclusionProbabilities)
			: calculateInclusionProbabilities(calculateInclusionProbabilities)
		{}
		//Should we calculate the inclusion probabilities? They're expensive. 
		bool calculateInclusionProbabilities;
		std::vector<mpfr_class> inclusionProbabilities;
	};
	void conditionalPoissonRejective(conditionalPoissonRejectiveArgs& args, boost::mt19937& randomSource);
}
#endif
