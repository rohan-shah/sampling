#ifndef SAMPLING_BASE_DOUBLE_PRECISION_HEADER_GUARD
#define SAMPLING_BASE_DOUBLE_PRECISION_HEADER_GUARD
#include <vector>
#include <algorithm>
namespace samplingDouble
{
	void sampfordBase(int n, std::vector<int>& indices, std::vector<double>& weights, std::vector<double>& rescaledWeights, std::vector<bool>& zeroWeights, std::vector<bool>& deterministicInclusion, int& nDeterministic, int& nZeroWeights);
	void samplingBase(int n, std::vector<int>& indices, std::vector<double>& weights, std::vector<bool>& zeroWeights, std::vector<bool>& deterministicInclusion, int& nDeterministic, int& nZeroWeights);
}
#endif
