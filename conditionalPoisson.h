#ifndef CONDITIONAL_POISSON_HEADER_GUARD
#define CONDITIONAL_POISSON_HEADER_GUARD
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
		{}
		std::size_t n;
		std::vector<bool> deterministicInclusion;
		std::vector<mpfr_class> exponentialParameters;
		std::vector<mpfr_class> expExponentialParameters;
		boost::numeric::ublas::matrix<mpfr_class> expNormalisingConstant;
		bool calculateInclusionProbabilities;
	};
	void conditionalPoisson(conditionalPoissonArgs& args, std::vector<int>& indices, std::vector<mpfr_class>& inclusionProbabilities, std::vector<mpfr_class>& weights, boost::mt19937& randomSource);
	void conditionalPoissonInclusionProbabilities(conditionalPoissonArgs& args, std::vector<mpfr_class>& inclusionProbabilities, std::vector<mpfr_class>& weights);
	void calculateExpNormalisingConstants(std::vector<mpfr_class>& expExponentialParameters, std::vector<mpfr_class>& exponentialParameters, boost::numeric::ublas::matrix<mpfr_class>& expNormalisingConstant, int n, int nUnits, std::vector<bool>& ignore);
}
#endif
