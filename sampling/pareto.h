#ifndef PARETO_SAMPLING_HEADER_GUARD
#define PARETO_SAMPLING_HEADER_GUARD
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include "includeMPFRSampling.h"
namespace sampling
{
	struct paretoSamplingArgs
	{
	public:
		paretoSamplingArgs ()
		{}
		//The indices of the sampled units
		std::vector<int> indices;
		//The inclusion probabilities. Contains nUnits entries, but only the ones corresponding to selected units are actually set.
		std::vector<mpfr_class> inclusionProbabilities;
		//A copy of the sampling weights made after deterministically selected units are removed, and then the sizes are rescaled.
		std::vector<mpfr_class> rescaledWeights;
		//The size variables to use for the sampling
		std::vector<mpfr_class> weights;
		//The number of units to select
		std::size_t n;
		std::vector<bool> deterministicInclusion;
		std::vector<bool> zeroWeights;
		struct paretoStatistic
		{
			double statistic;
			int order;
			bool operator<(const paretoStatistic& other)
			{
				return statistic < other.statistic;
			}
		};
		std::vector<paretoStatistic> paretoStatistics;
	};
	void pareto(paretoSamplingArgs& args, boost::mt19937& randomSource);
}
#endif
