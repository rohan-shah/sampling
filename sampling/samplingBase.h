#ifndef SAMPLING_BASE_HEADER_GUARD
#define SAMPLING_BASE_HEADER_GUARD
#include <vector>
#include <algorithm>
#include "includeMPFRSampling.h"
namespace sampling
{
	void sampfordBase(int n, std::vector<int>& indices, std::vector<mpfr_class>& weights, std::vector<mpfr_class>& rescaledWeights, std::vector<bool>& zeroWeights, std::vector<bool>& deterministicInclusion, int& nDeterministic, int& nZeroWeights);
	void samplingBase(int n, std::vector<int>& indices, std::vector<mpfr_class>& weights, std::vector<bool>& zeroWeights, std::vector<bool>& deterministicInclusion, int& nDeterministic, int& nZeroWeights);
}
#endif
