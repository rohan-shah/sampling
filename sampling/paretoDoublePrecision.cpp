#include "paretoDoublePrecision.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "samplingBaseDoublePrecision.h"
#include <functional>
namespace samplingDouble
{
	void pareto(paretoSamplingArgs& args, boost::mt19937& randomSource)
	{
		std::vector<int>& indices = args.indices;
		std::vector<double>& weights = args.weights;
		int nUnits = (int)weights.size();

		int nZeroWeights = 0, nDeterministic = 0;
		samplingBase(args.n, indices, weights, args.zeroWeights, args.deterministicInclusion, nDeterministic, nZeroWeights);
		if(nDeterministic == (int)args.n) return;

		boost::random::uniform_real_distribution<> standardUniform(0, 1);
		//Now compute the pareto statistics
		args.paretoStatistics.clear();
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i] && !args.zeroWeights[i])
			{
				double uniform = standardUniform(randomSource);
				paretoSamplingArgs::paretoStatistic newStatistic;
				newStatistic.statistic = (uniform * (1 - weights[i]))/(weights[i]*(1-uniform));
				newStatistic.order = i;
				args.paretoStatistics.push_back(newStatistic);
			}
		}
		std::sort(args.paretoStatistics.begin(), args.paretoStatistics.end());
		//Select so many smallest values
		for(int i = 0; i < (int)(args.n - nDeterministic); i++)
		{
			indices.push_back(args.paretoStatistics[i].order);
		}
	}
}
