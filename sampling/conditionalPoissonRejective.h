#ifndef CONDITIONAL_POISSON_REJECTIVE_HEADER_GUARD
#define CONDITIONAL_POISSON_REJECTIVE_HEADER_GUARD
#include <vector>
#include "includeMPFRSampling.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/numeric/ublas/matrix.hpp>
namespace sampling
{
	struct conditionalPoissonArgs
	{
	public:
		conditionalPoissonArgs()
			: indices(NULL), inclusionProbabilities(NULL), weights(NULL), rescaledWeights(NULL)
		{}
		//The vector of selected units
		std::vector<int>* indices;
		//The inclusion probabilities. Contains nUnits entries, but only the ones corresponding to selected units are actually set.
		std::vector<mpfr_class>* inclusionProbabilities;
		//The input size variables for the sampling
		std::vector<mpfr_class>* weights;
		//A copy of the sampling weights made after deterministically selected units are removed, and then the sizes are rescaled.
		std::vector<mpfr_class>* rescaledWeights;
		//Number of units to select
		std::size_t n;
		//Vector indicating which units were deterministically selected
		std::vector<bool> deterministicInclusion;
		//Vector indicating which units had zero weight
		std::vector<bool> zeroWeights;
		//Should we calculate the inclusion probabilities? They're expensive. 
		bool calculateInclusionProbabilities;

		//All the remaining members are internal working memory
		
		//Units not subject to the usual sampling procedure. Essentially deterministicInclusion & zeroWeights
		std::vector<bool> ignore;
		std::vector<mpfr_class> exponentialParameters;
		std::vector<mpfr_class> expExponentialParameters;
		boost::numeric::ublas::matrix<mpfr_class> expNormalisingConstant;
	private:
		conditionalPoissonArgs(const conditionalPoissonArgs& other);
		conditionalPoissonArgs& operator=(const conditionalPoissonArgs& other);
	};
	void conditionalPoissonRejective(conditionalPoissonArgs& args, boost::mt19937& randomSource);
	void conditionalPoissonInclusionProbabilities(conditionalPoissonArgs& args);
	void calculateExpNormalisingConstants(std::vector<mpfr_class>& expExponentialParameters, std::vector<mpfr_class>& exponentialParameters, boost::numeric::ublas::matrix<mpfr_class>& expNormalisingConstant, int n, int nUnits, std::vector<bool>& ignore);
}
#endif
