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
	void sampfordFromParetoNaive(sampfordFromParetoNaiveArgs& args, std::vector<int>& indices, std::vector<mpfr_class>& inclusionProbabilities, const std::vector<mpfr_class>& weights, boost::mt19937& randomSource, std::vector<mpfr_class>& copiedWeights)
	{
		int nUnits = weights.size();
		args.paretoArgs.n = args.n;
		args.paretoArgs.calculateInclusionProbabilities = false;
		pareto(args.paretoArgs, indices, inclusionProbabilities, weights, randomSource, copiedWeights);
		inclusionProbabilities.resize(nUnits);
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.paretoArgs.deterministicInclusion[i])
			{
				inclusionProbabilities[i] = copiedWeights[i];
			}
			else inclusionProbabilities[i] = 1;
		}
	}
}
