#include "conditionalPoissonRejective.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "samplingBase.h"
namespace sampling
{
	using boost::multiprecision::log;
	using std::log;
	using boost::multiprecision::exp;
	using std::exp;
	void conditionalPoissonRejective(conditionalPoissonArgs& args, boost::mt19937& randomSource)
	{
		std::vector<int>& indices = *args.indices;
		std::vector<mpfr_class>& weights = *args.weights;
		std::vector<mpfr_class>& inclusionProbabilities = *args.inclusionProbabilities;
		std::vector<mpfr_class>& rescaledWeights = *args.rescaledWeights;
		int nUnits = (int)weights.size();

		if(samplingBase(args.n, indices, weights, rescaledWeights, args.zeroWeights, args.deterministicInclusion, inclusionProbabilities)) return;
		int deterministicIndices = indices.size();
		computeExponentialParameters(args);
beginSample:
		indices.resize(deterministicIndices);
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i])
			{
				boost::random::bernoulli_distribution<double> currentUnitDist((double)rescaledWeights[i]);
				if(currentUnitDist(randomSource))
				{
					indices.push_back(i);
				}
			}
			//If we can't reach the target of n units, then start again. 
			if(indices.size() > args.n || indices.size() + (nUnits - i - 1) < args.n) goto beginSample;
		}
		if(indices.size() != args.n)
		{
			throw std::runtime_error("Internal error");
		}
		if(args.calculateInclusionProbabilities) conditionalPoissonInclusionProbabilities(args);
	}
}
