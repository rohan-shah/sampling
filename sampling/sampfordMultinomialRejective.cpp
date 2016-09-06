#include "sampford.h"
#include <boost/random/uniform_real_distribution.hpp>
#include <functional>
#include "samplingBase.h"
namespace sampling
{
	void sampfordMultinomialRejective(sampfordMultinomialRejectiveArgs& args, boost::mt19937& randomSource)
	{
		std::vector<int>& indices = args.indices;
		std::vector<mpfr_class>& weights = args.weights;
		std::vector<mpfr_class>& rescaledWeights = args.rescaledWeights;
		std::vector<bool>& zeroWeights = args.zeroWeights;
		std::vector<bool>& deterministicInclusion = args.deterministicInclusion;
		std::vector<int> deterministicIndices;
		int nDeterministic = 0, nZeroWeights = 0;
		sampfordBase((int)args.n, deterministicIndices, weights, rescaledWeights, zeroWeights, deterministicInclusion, nDeterministic, nZeroWeights);
		int nUnits = (int)weights.size();
		if((int)args.n == nDeterministic)
		{
			indices.swap(deterministicIndices);
			return;
		}
		args.accumulated1.clear();
		args.accumulated2.clear();
		args.accumulated1.reserve(nUnits);
		args.accumulated2.reserve(nUnits);
		args.indices1.clear();

		//Rescale weights so they sum to args.n. These are the weights used for the first sample
		mpfr_class rescaledCumulative = 0;
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i] && !args.zeroWeights[i])
			{
				mpfr_class& prob = rescaledWeights[i];
				rescaledCumulative += prob;
				args.accumulated1.push_back((double)rescaledCumulative);
				args.indices1.push_back(i);
			}
		}

		//Accumulated weights for the other samples
		mpfr_class cumulativeOther = 0;
		for(int i = 0; i < nUnits; i++)
		{
			if(!zeroWeights[i] && !deterministicInclusion[i])
			{
				mpfr_class& rescaled = rescaledWeights[i];
				mpfr_class asProb = rescaled / (1 - rescaled);
				cumulativeOther += asProb;
				args.accumulated2.push_back((double)cumulativeOther);
			}
		}
		boost::random::uniform_real_distribution<> firstSampleDist(0, (double)rescaledCumulative);
getSample:
		indices = deterministicIndices;
		double firstSample = firstSampleDist(randomSource);
		int firstIndex = (int)std::distance(args.accumulated1.begin(), std::upper_bound(args.accumulated1.begin(), args.accumulated1.end(), firstSample, std::less_equal<double>()));
		indices.push_back(args.indices1[firstIndex]);
		{
			boost::random::uniform_real_distribution<> otherSamplesDist(0, (double)cumulativeOther);
			for(int i = 0; i < (int)args.n-nDeterministic-1; i++)
			{
				double sampledWeight = otherSamplesDist(randomSource);
				int index = (int)std::distance(args.accumulated2.begin(), std::upper_bound(args.accumulated2.begin(), args.accumulated2.end(), sampledWeight, std::less_equal<double>()));
				//Could in theory happen if we get exactly the right value from the RNG
				if(args.deterministicInclusion[index] || args.zeroWeights[index])
				{
					goto getSample;
				}
				indices.push_back(args.indices1[index]);
			}
			std::sort(indices.begin(), indices.end());
			if(std::unique(indices.begin(), indices.end()) != indices.end())
			{
				goto getSample;
			}
		}
	}
}
