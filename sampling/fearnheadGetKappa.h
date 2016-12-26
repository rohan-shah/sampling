#ifndef FEARNHEAD_GET_KAPPA_HEADER_GUARD
#define FEARNHEAD_GET_KAPPA_HEADER_GUARD
#include "includeMPFRSampling.h"
#include <vector>
#include <boost/random/mersenne_twister.hpp>
namespace samplingDouble
{
	void fearnheadGetKappa(std::vector<double>& sortedWeights, boost::mt19937& randomSource, int N, int& A, double& B);
}
#endif

