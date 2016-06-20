#include "conditionalPoissonRejective.h"
#include <boost/random/bernoulli_distribution.hpp>
#include "samplingBase.h"
namespace sampling
{
	using boost::multiprecision::log;
	using std::log;
	using boost::multiprecision::exp;
	using std::exp;
	void conditionalPoissonRejective(conditionalPoissonRejectiveArgs& args, boost::mt19937& randomSource)
	{
		std::vector<int>& indices = args.indices;
		indices.clear();
		std::vector<mpfr_class>& weights = args.weights;
		std::vector<mpfr_class>& inclusionProbabilities = args.inclusionProbabilities;
		std::vector<bool>& deterministicInclusion = args.deterministicInclusion;
		std::vector<bool>& zeroWeights = args.zeroWeights;
		int nUnits = (int)weights.size();
		int nZeroWeights = 0, nDeterministic = 0;
		samplingBase(args.n, indices, weights, zeroWeights, deterministicInclusion, nDeterministic, nZeroWeights);
		computeExponentialParameters(args);
beginSample:
		int copiedZeroWeights = nZeroWeights, copiedDeterministic = nDeterministic;
		indices.resize(nDeterministic);
		for(int i = 0; i < nUnits; i++)
		{
			if(zeroWeights[i]) copiedZeroWeights--;
			else if(deterministicInclusion[i]) copiedDeterministic--;
			else
			{
				boost::random::bernoulli_distribution<double> currentUnitDist((double)weights[i]);
				if(currentUnitDist(randomSource))
				{
					indices.push_back(i);
				}
				//If we can't reach the target of n units, then start again. 
				if(indices.size() > args.n || indices.size() + (nUnits - i - 1) - copiedZeroWeights - copiedDeterministic < args.n) goto beginSample;
			}
		}
		if(indices.size() != args.n)
		{
			throw std::runtime_error("Internal error");
		}
		if(args.calculateInclusionProbabilities) conditionalPoissonInclusionProbabilities(args, inclusionProbabilities);
	}
}
