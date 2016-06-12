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
		std::vector<int>& indices = *args.indices;
		std::vector<mpfr_class>& weights = *args.weights;
		std::vector<mpfr_class>& rescaledWeights = *args.rescaledWeights;
		std::vector<mpfr_class>& inclusionProbabilities = *args.inclusionProbabilities;
		int nUnits = (int)weights.size();
		indices.clear();

		if(samplingBase(args.n, indices, weights, rescaledWeights, args.zeroWeights, args.deterministicInclusion, inclusionProbabilities)) return;
		int deterministicIndices = indices.size();
		computeExponentialParameters(args);
		conditionalPoissonInclusionProbabilities(args, inclusionProbabilities);
		int chosen = 0;
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i] && !args.zeroWeights[i])
			{
				double parameter;
				if(args.n - deterministicIndices - 1 - chosen == 0)
				{
					parameter = (args.expExponentialParameters[i] / args.expNormalisingConstant(i, args.n - deterministicIndices - chosen - 1)).convert_to<double>();
				}
				else
				{
					parameter = (args.expExponentialParameters[i] * args.expNormalisingConstant(i+1, args.n - deterministicIndices - 1 - chosen - 1) / args.expNormalisingConstant(i, args.n - deterministicIndices - chosen - 1)).convert_to<double>();
				}
				boost::random::bernoulli_distribution<> bernoulli(parameter);
				if(bernoulli(randomSource))
				{
					indices.push_back(i);
					chosen++;
				}
			}
			if(chosen == (int)args.n - deterministicIndices) break;
		}
	}
}
