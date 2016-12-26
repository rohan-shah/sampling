#ifndef FEARNHEAD_SAMPLING_HEADER_GUARD
#define FEARNHEAD_SAMPLING_HEADER_GUARD
#include "includeMPFRSampling.h"
#include <vector>
#include <boost/random/mersenne_twister.hpp>
namespace samplingDouble
{
	struct fearnheadSamplingArgs
	{
		int n;
		double c;
		std::vector<double> weights;
		std::vector<double> sortedWeights;
		std::vector<int> indices;
		std::vector<bool> deterministicInclusion;
	};
	void fearnheadSampling(fearnheadSamplingArgs& args, boost::mt19937& randomSource);
}
namespace sampling
{
	struct fearnheadSamplingArgs
	{
		int n;
		mpfr_class c;
		std::vector<mpfr_class> weights;
		std::vector<mpfr_class> sortedWeights;
		std::vector<int> indices;
		std::vector<bool> deterministicInclusion;
	};
	void fearnheadSampling(fearnheadSamplingArgs& args, boost::mt19937& randomSource);
}
#endif
