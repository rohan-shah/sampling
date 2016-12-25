#ifndef SYSTEMATIC_SAMPLING_HEADER_GUARD
#define SYSTEMATIC_SAMPLING_HEADER_GUARD
#include "includeMPFRSampling.h"
#include <vector>
#include <boost/random/mersenne_twister.hpp>
namespace sampling
{
	void systematicSamplingDouble(const std::vector<double>& weights, double interval, std::vector<int>& indices, boost::mt19937& randomSource);
	void systematicSampling(const std::vector<mpfr_class>& weights, mpfr_class interval, std::vector<int>& indices, boost::mt19937& randomSource);
}
#endif
