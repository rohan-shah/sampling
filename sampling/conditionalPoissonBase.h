#ifndef CONDITIONAL_POISSON_BASE_HEADER_GUARD
#define CONDITIONAL_POISSON_BASE_HEADER_GUARD
#include <vector>
#include <algorithm>
#include <boost/numeric/ublas/matrix.hpp>
#include "includeMPFRSampling.h"
namespace sampling
{
	struct conditionalPoissonArgs
	{
	public:
		conditionalPoissonArgs()
		{}
		//The vector of selected units
		std::vector<int> indices;
		std::vector<mpfr_class> weights;
		//Number of units to select
		std::size_t n;

		std::vector<mpfr_class> exponentialParameters;
		std::vector<mpfr_class> expExponentialParameters;
		boost::numeric::ublas::matrix<mpfr_class> expNormalisingConstant;
		std::vector<bool> deterministicInclusion, zeroWeights;
	private:
		conditionalPoissonArgs(const conditionalPoissonArgs& other);
		conditionalPoissonArgs& operator=(const conditionalPoissonArgs& other);
	};
	void conditionalPoissonInclusionProbabilities(conditionalPoissonArgs& args, std::vector<mpfr_class>& inclusionProbabilities);
	void conditionalPoissonSecondOrderInclusionProbabilities(conditionalPoissonArgs& args, std::vector<mpfr_class>& inclusionProbabilities, boost::numeric::ublas::matrix<mpfr_class>& secondOrder);
	void calculateExpNormalisingConstants(std::vector<mpfr_class>& expExponentialParameters, std::vector<mpfr_class>& exponentialParameters, boost::numeric::ublas::matrix<mpfr_class>& expNormalisingConstant, int n, int nUnits, std::vector<bool>& zeroWeights, std::vector<bool>& deterministicInclusion);
	void computeExponentialParameters(conditionalPoissonArgs& args);
	void calculateExpNormalisingConstants(conditionalPoissonArgs& args);
}
#endif
