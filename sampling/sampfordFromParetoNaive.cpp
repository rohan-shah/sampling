#include "sampford.h"
#include "pareto.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <functional>
namespace sampling
{
	using boost::multiprecision::log;
	using std::log;
	using boost::multiprecision::exp;
	using std::exp;
	void sampfordFromParetoNaive(sampfordFromParetoNaiveArgs& args, boost::mt19937& randomSource)
	{
		args.paretoArgs.n = args.n;
		std::swap(args.paretoArgs.weights, args.weights);

		pareto(args.paretoArgs, randomSource);

		std::swap(args.paretoArgs.weights, args.weights);
		std::swap(args.paretoArgs.rescaledWeights, args.rescaledWeights);
		std::swap(args.paretoArgs.indices, args.indices);
		std::swap(args.paretoArgs.deterministicInclusion, args.deterministicInclusion);
		std::swap(args.paretoArgs.zeroWeights, args.zeroWeights);

		std::vector<mpfr_class>& weights = args.weights;
		std::vector<mpfr_class>& inclusionProbabilities = args.inclusionProbabilities;
		std::vector<mpfr_class>& rescaledWeights = args.rescaledWeights;
		int nUnits = weights.size();
		inclusionProbabilities.resize(nUnits);
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i])
			{
				inclusionProbabilities[i] = rescaledWeights[i];
			}
			else 
			{
				if(weights[i] == 0) inclusionProbabilities[i] = 0;
				else inclusionProbabilities[i] = 1;
			}
		}
	}
}
