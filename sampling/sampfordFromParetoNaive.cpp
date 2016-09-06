#include "sampford.h"
#include "pareto.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <functional>
#include "samplingBase.h"
namespace sampling
{
	using boost::multiprecision::log;
	using std::log;
	using boost::multiprecision::exp;
	using std::exp;
	void sampfordFromParetoNaive(sampfordFromParetoNaiveArgs& args, boost::mt19937& randomSource)
	{
		std::vector<mpfr_class>& weights = args.weights;
		std::vector<mpfr_class>& rescaledWeights = args.rescaledWeights;
		int nDeterministic = 0, nZeroWeights = 0;
		sampfordBase((int)args.n, args.indices, weights, rescaledWeights, args.zeroWeights, args.deterministicInclusion, nDeterministic, nZeroWeights);
		if(args.indices.size() == args.n)
		{
			return;
		}
		args.paretoArgs.n = args.n;
		std::swap(args.paretoArgs.weights, args.rescaledWeights);

		pareto(args.paretoArgs, randomSource);

		std::swap(args.paretoArgs.weights, args.rescaledWeights);
		std::swap(args.paretoArgs.indices, args.indices);
		std::swap(args.paretoArgs.deterministicInclusion, args.deterministicInclusion);
		std::swap(args.paretoArgs.zeroWeights, args.zeroWeights);
	}
}
