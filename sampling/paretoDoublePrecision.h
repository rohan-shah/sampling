#ifndef PARETO_SAMPLING_DOUBLE_PRECISION_HEADER_GUARD
#define PARETO_SAMPLING_DOUBLE_PRECISION_HEADER_GUARD
#include <vector>
#include <boost/random/mersenne_twister.hpp>
namespace samplingDouble
{
	struct paretoSamplingArgs
	{
	public:
		paretoSamplingArgs ()
		{}
		//The indices of the sampled units
		std::vector<int> indices;
		//The size variables to use for the sampling
		std::vector<double> weights;
		//The number of units to select
		std::size_t n;
		std::vector<bool> deterministicInclusion;
		std::vector<bool> zeroWeights;
		struct paretoStatistic
		{
			double statistic;
			int order;
			bool operator<(const paretoStatistic& other) const
			{
				return statistic < other.statistic;
			}
		};
		std::vector<paretoStatistic> paretoStatistics;
	};
	void pareto(paretoSamplingArgs& args, boost::mt19937& randomSource);
}
#endif
