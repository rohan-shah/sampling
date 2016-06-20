#include "sampford.h"
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
	void sampfordConditionalPoissonRejective(sampfordConditionalPoissonRejectiveArgs& args, boost::mt19937& randomSource)
	{
		std::vector<int>& indices = args.indices;
		std::vector<mpfr_class>& weights = args.weights;
		int nUnits = (int)weights.size();
		std::vector<mpfr_class>& rescaledWeights = args.rescaledWeights;
		std::vector<bool>& deterministicInclusion = args.deterministicInclusion;
		std::vector<bool>& zeroWeights = args.zeroWeights;
		int nDeterministic = 0, nZeroWeights = 0;
		sampfordBase(args.n, indices, weights, rescaledWeights, zeroWeights, deterministicInclusion, nDeterministic, nZeroWeights);
		if((int)args.n == nDeterministic)
		{
			return;
		}
		args.accumulated.clear();
		args.accumulated.reserve(nUnits - nDeterministic);
		args.indices1.clear();
		mpfr_class cumulative = 0;
		for(int i = 0; i < nUnits; i++)
		{
			if(!zeroWeights[i] && !deterministicInclusion[i])
			{
				cumulative += weights[i];
				args.accumulated.push_back((double)cumulative);
				args.indices1.push_back(i);
			}
		}
		boost::random::uniform_real_distribution<> firstSampleDist(0, (double)cumulative);
beginSample:
		indices.resize(nDeterministic);
		double firstSample = firstSampleDist(randomSource);
		int firstIndex = (int)std::distance(args.accumulated.begin(), std::upper_bound(args.accumulated.begin(), args.accumulated.end(), firstSample, std::less_equal<double>()));
		firstIndex = args.indices1[firstIndex];
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i] && !args.zeroWeights[i])
			{
				boost::random::bernoulli_distribution<double> currentUnitDist((double)rescaledWeights[i]);
				if(currentUnitDist(randomSource))
				{
					indices.push_back(i);
				}
			}
			//If we can't reach the target of n units, then start again. 
			if(indices.size() > args.n-1 || indices.size() + (nUnits - i - 1) < args.n-1) goto beginSample;
		}
		if(indices.size() != args.n-1)
		{
			throw std::runtime_error("Internal error");
		}
		if(args.n > 1 && std::find(indices.begin(), indices.end(), firstIndex) != indices.end())
		{
			goto beginSample;
		}
		indices.push_back(firstIndex);
		std::sort(indices.begin(), indices.end());
		if(std::unique(indices.begin(), indices.end()) != indices.end())
		{
			throw std::runtime_error("Internal error");
		}
	}
}
