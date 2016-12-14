#ifndef SAMPFORD_SAMPLING_DOUBLE_PRECISION_HEADER_GUARD
#define SAMPFORD_SAMPLING_DOUBLE_PRECISION_HEADER_GUARD
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include "paretoDoublePrecision.h"
namespace samplingDouble
{
	struct sampfordFromParetoNaiveArgs
	{
	public:
		paretoSamplingArgs paretoArgs;
		sampfordFromParetoNaiveArgs()
		{}
		std::size_t n;
		std::vector<double> weights;
		std::vector<double> rescaledWeights;
		std::vector<int> indices;
		std::vector<bool> deterministicInclusion;
		std::vector<bool> zeroWeights;
	};
	void sampfordFromParetoNaive(sampfordFromParetoNaiveArgs& args, boost::mt19937& randomSource);

}
#endif
