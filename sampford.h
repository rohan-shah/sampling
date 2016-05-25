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
		std::size_t n;
		std::vector<double> accumulated1;
		std::vector<double> accumulated2;
		std::vector<bool> deterministicInclusion;
		std::vector<int> deterministicIndices;
	};
	void sampfordMultinomialRejective(sampfordMultinomialRejectiveArgs& args, std::vector<int>& indices, std::vector<mpfr_class>& inclusionProbabilities, std::vector<mpfr_class>& weights, boost::mt19937& randomSource);
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
	};
	void sampfordFromParetoNaive(sampfordFromParetoNaiveArgs& args, std::vector<int>& indices, std::vector<mpfr_class>& inclusionProbabilities, std::vector<mpfr_class>& weights, boost::mt19937& randomSource);

}
#endif
