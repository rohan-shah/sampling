#include "conditionalPoissonSequential.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "samplingBase.h"
namespace sampling
{
	using boost::multiprecision::log;
	using std::log;
	using boost::multiprecision::exp;
	using std::exp;
	void conditionalPoissonSequential(conditionalPoissonSequentialArgs& args, boost::mt19937& randomSource)
	{
		std::vector<int>& indices = args.indices;
		std::vector<mpfr_class>& weights = args.weights;
		std::vector<mpfr_class>& inclusionProbabilities = args.inclusionProbabilities;
		std::vector<bool>& zeroWeights = args.zeroWeights;
		std::vector<bool>& deterministicInclusion = args.deterministicInclusion;
		int nUnits = (int)weights.size();
		indices.clear();

		int nZeroWeights = 0, nDeterministic = 0;
		samplingBase(args.n, indices, weights, zeroWeights, deterministicInclusion, nDeterministic, nZeroWeights);

		computeExponentialParameters(args);
		if(args.calculateInclusionProbabilities)
		{
			conditionalPoissonInclusionProbabilities(args, inclusionProbabilities);
		}
		else calculateExpNormalisingConstants(args);
		int chosen = 0;
		int skipped = 0;
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.zeroWeights[i] && !args.deterministicInclusion[i])
			{
				double parameter;
				if(args.n - nDeterministic - 1 - chosen == 0)
				{
					parameter = mpfr_class(args.expExponentialParameters[i] / args.expNormalisingConstant(i-skipped, args.n - nDeterministic - chosen - 1)).convert_to<double>();
				}
				else
				{
					parameter = mpfr_class(args.expExponentialParameters[i] * args.expNormalisingConstant(i+1-skipped, args.n - nDeterministic - 1 - chosen - 1) / args.expNormalisingConstant(i-skipped, args.n - nDeterministic - chosen - 1)).convert_to<double>();
				}
#ifndef NDEBUG
				if(parameter > 1) throw std::runtime_error("Internal error");
#endif
				boost::random::bernoulli_distribution<> bernoulli(parameter);
				if(bernoulli(randomSource))
				{
					indices.push_back(i);
					chosen++;
				}
			}
			else skipped++;
			if(chosen == (int)args.n - nDeterministic) break;
		}
	}
}
