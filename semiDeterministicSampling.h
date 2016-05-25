#ifndef SEMIDETERMINISTIC_SAMPLING_HEADER_GUARD
#define SEMIDETERMINISTIC_SAMPLING_HEADER_GUARD
#include <vector>
#include "includeMPFRSampling.h"
#include <boost/random/mersenne_twister.hpp>
namespace sampling
{
	struct semiDeterministicSamplingArgs
	{
	public:
		semiDeterministicSamplingArgs()
		{}
		std::size_t n;
		std::vector<bool> deterministicInclusion;
		std::vector<int> deterministicIndices;
	};
	void semiDeterministic(semiDeterministicSamplingArgs& args, std::vector<int>& indices, std::vector<mpfr_class>& inclusionProbabilities, std::vector<mpfr_class>& weights, boost::mt19937& randomSource);
}
#endif
