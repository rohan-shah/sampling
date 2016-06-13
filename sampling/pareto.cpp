#include "pareto.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "samplingBase.h"
#include <functional>
namespace sampling
{
	using boost::multiprecision::log;
	using std::log;
	using boost::multiprecision::exp;
	using std::exp;
	void pareto(paretoSamplingArgs& args, boost::mt19937& randomSource)
	{
		std::vector<int>& indices = args.indices;
		std::vector<mpfr_class>& weights = args.weights;
		std::vector<mpfr_class>& inclusionProbabilities = args.inclusionProbabilities;
		std::vector<mpfr_class>& rescaledWeights = args.rescaledWeights;
		int nUnits = (int)weights.size();

		if(samplingBase(args.n, indices, weights, rescaledWeights, args.zeroWeights, args.deterministicInclusion, inclusionProbabilities)) return;
		int deterministicIndices = indices.size();

		boost::random::uniform_real_distribution<> standardUniform(0, 1);
		//Now compute the pareto statistics
		args.paretoStatistics.clear();
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i])
			{
				double uniform = standardUniform(randomSource);
				paretoSamplingArgs::paretoStatistic newStatistic;
				mpfr_class value = (uniform * (1 - rescaledWeights[i]))/(rescaledWeights[i]*(1-uniform));
				newStatistic.statistic = value.convert_to<double>();
				newStatistic.order = i;
				args.paretoStatistics.push_back(newStatistic);
			}
		}
		std::sort(args.paretoStatistics.begin(), args.paretoStatistics.end());
		//Select so many smallest values
		for(int i = 0; i < (int)(args.n - deterministicIndices); i++)
		{
			indices.push_back(args.paretoStatistics[i].order);
		}
	}
}
