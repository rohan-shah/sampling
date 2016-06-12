#ifndef SAMPLING_BASE_HEADER_GUARD
#define SAMPLING_BASE_HEADER_GUARD
#include <vector>
#include <algorithm>
#include "includeMPFRSampling.h"
namespace sampling
{
	bool samplingBase(int n, std::vector<int>& indices, std::vector<mpfr_class>& weights, std::vector<mpfr_class>& rescaledWeights, std::vector<bool>& zeroWeights, std::vector<bool>& deterministicInclusion, std::vector<mpfr_class>& inclusionProbabilities);
}
#endif
