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
		std::vector<mpfr_class>& weights = *args.weights;
		int nUnits = weights.size();
		args.paretoArgs.n = args.n;
		args.paretoArgs.calculateInclusionProbabilities = false;
		args.paretoArgs.inclusionProbabilities = args.inclusionProbabilities;
		args.paretoArgs.weights = args.weights;
		args.paretoArgs.rescaledWeights = args.rescaledWeights;
		args.paretoArgs.indices = args.indices;

		pareto(args.paretoArgs, randomSource);

		std::vector<mpfr_class>& inclusionProbabilities = *args.inclusionProbabilities;
		std::vector<mpfr_class>& rescaledWeights= *args.rescaledWeights;
		inclusionProbabilities.resize(nUnits);
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.paretoArgs.deterministicInclusion[i])
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
