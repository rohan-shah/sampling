#ifndef SAMPFORD_SAMPLING_HEADER_GUARD
#define SAMPFORD_SAMPLING_HEADER_GUARD
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include "includeMPFRSampling.h"
#include "pareto.h"
namespace sampling
{
	struct sampfordMultinomialRejectiveArgs
	{
	public:
		sampfordMultinomialRejectiveArgs()
		{}
		//A copy of the sampling weights made after deterministically selected units are removed, and then the sizes are rescaled.
		std::vector<mpfr_class> rescaledWeights;
		//The size variables for the sampling
		std::vector<mpfr_class> weights;
		//Inclusion probabilities. Contains nUnits values, but only the values for the chosen units are set. 
		std::vector<mpfr_class> inclusionProbabilities;
		//Indices of selected units
		std::vector<int> indices;
		//Number of units to sample
		std::size_t n;
		//Working memory
		std::vector<double> accumulated1;
		std::vector<double> accumulated2;
		std::vector<bool> deterministicInclusion;
		std::vector<int> deterministicIndices;
	};
	void sampfordMultinomialRejective(sampfordMultinomialRejectiveArgs& args, boost::mt19937& randomSource);
	struct sampfordConditionalPoissonRejectiveArgs
	{
	public:
		sampfordConditionalPoissonRejectiveArgs()
		{}
		std::size_t n;
		std::vector<bool> deterministicInclusion;
		std::vector<double> accumulated;
	};
	void sampfordConditionalPoissonRejective(sampfordConditionalPoissonRejectiveArgs& args, std::vector<int>& indices, std::vector<mpfr_class>& inclusionProbabilities, std::vector<mpfr_class>& weights, boost::mt19937& randomSource);
	struct sampfordFromParetoNaiveArgs
	{
	public:
		paretoSamplingArgs paretoArgs;
		sampfordFromParetoNaiveArgs()
		{}
		std::size_t n;
		std::vector<mpfr_class> weights;
		std::vector<mpfr_class> rescaledWeights;
		std::vector<mpfr_class> inclusionProbabilities;
		std::vector<int> indices;
		std::vector<bool> deterministicInclusion;
		std::vector<bool> zeroWeights;
	};
	void sampfordFromParetoNaive(sampfordFromParetoNaiveArgs& args, boost::mt19937& randomSource);

}
#endif
