#include "conditionalPoissonBase.h"
namespace sampling
{
	void calculateExpNormalisingConstants(conditionalPoissonArgs& args)
	{
		std::vector<mpfr_class>& rescaledWeights = args.rescaledWeights;
		std::vector<bool>& ignore = args.ignore;
		int nUnits = (int)rescaledWeights.size();
		int deterministicIndices = 0, ignoreIndices = 0;
		ignore.resize(nUnits);
		for(int i = 0; i < nUnits; i++)
		{
			ignore[i] = args.deterministicInclusion[i] || args.zeroWeights[i];
			if(args.deterministicInclusion[i]) deterministicIndices++;
			if(ignore[i]) ignoreIndices++;
		}
		//Now compute the inclusion probabilities
		calculateExpNormalisingConstants(args.expExponentialParameters, args.exponentialParameters, args.expNormalisingConstant, (int)args.n - deterministicIndices, nUnits - ignoreIndices, args.ignore);
	}
	void conditionalPoissonInclusionProbabilities(conditionalPoissonArgs& args, std::vector<mpfr_class>& inclusionProbabilities)
	{
		std::vector<mpfr_class>& rescaledWeights = args.rescaledWeights;
		std::vector<bool>& ignore = args.ignore;
		int nUnits = (int)rescaledWeights.size();
		int deterministicIndices = 0, ignoreIndices = 0;
		ignore.resize(nUnits);
		for(int i = 0; i < nUnits; i++)
		{
			ignore[i] = args.deterministicInclusion[i] || args.zeroWeights[i];
			if(args.deterministicInclusion[i]) deterministicIndices++;
			if(ignore[i]) ignoreIndices++;
		}
		//Now compute the inclusion probabilities
		inclusionProbabilities.resize(nUnits);
		calculateExpNormalisingConstants(args.expExponentialParameters, args.exponentialParameters, args.expNormalisingConstant, (int)args.n - deterministicIndices, nUnits - ignoreIndices, args.ignore);
		if((int)args.n == deterministicIndices) return;
		mpfr_class expNormalisingConstant = args.expNormalisingConstant(0, args.n - deterministicIndices - 1);
		for(int unitCounter = 0; unitCounter < nUnits; unitCounter++)
		{
			if(args.deterministicInclusion[unitCounter] || args.zeroWeights[unitCounter]) continue;
			//First term of the alternating sum
			if((args.n-deterministicIndices)% 2 == 1)
			{
				inclusionProbabilities[unitCounter] = 1;
			}
			else
			{
				inclusionProbabilities[unitCounter] = -1;
			}
			for(int j = 2; j <= (int)args.n-deterministicIndices; j++)
			{
				if((args.n - deterministicIndices- j) % 2 == 1)
				{
					inclusionProbabilities[unitCounter] -= args.expNormalisingConstant(0, j-2) / exp((j - 1)*args.exponentialParameters[unitCounter]);
				}
				else
				{
					inclusionProbabilities[unitCounter] += args.expNormalisingConstant(0, j-2) / exp((j - 1)*args.exponentialParameters[unitCounter]);
				}
			}
			inclusionProbabilities[unitCounter] *= exp((args.n - deterministicIndices) * args.exponentialParameters[unitCounter]) / expNormalisingConstant;
			if(inclusionProbabilities[unitCounter] <= 0)
			{
				std::stringstream ss;
				ss << "Inclusion probability had negative value " << inclusionProbabilities[unitCounter] << ", probably because of numerical instability";
				throw std::runtime_error(ss.str().c_str());
			}
		}
	}
	void calculateExpNormalisingConstants(std::vector<mpfr_class>& expExponentialParameters, std::vector<mpfr_class>& exponentialParameters, boost::numeric::ublas::matrix<mpfr_class>& expNormalisingConstant, int n, int nUnits, std::vector<bool>& ignore)
	{
		//We start by computing the normalising constants. First index is k, second is z. All indices in this loop are 1 indexed
		expNormalisingConstant.resize(nUnits,n);
		//This will skip over the *ignored* units (the ones that were deterministically included)
		int k = (int)expExponentialParameters.size();
		for(int unitIndex = nUnits; unitIndex >= 1; unitIndex--)
		{
			while(ignore[k-1]) k--;
			for(int z = 1; z <= std::min(n, nUnits-unitIndex+1); z++)
			{
				if(z == 1)
				{
					mpfr_class sum = 0;
					for(int unitIndex2 = k; unitIndex2 <= (int)expExponentialParameters.size(); unitIndex2++)
					{
						if(!ignore[unitIndex2-1]) sum += expExponentialParameters[unitIndex2-1];
					}
					expNormalisingConstant(unitIndex-1, z-1) = sum;
#ifndef NDEBUG
					assert(expNormalisingConstant(unitIndex-1, z-1) == expNormalisingConstant(unitIndex-1, z-1));
#endif
				}
				else if(z == nUnits - unitIndex + 1)
				{
					mpfr_class sum = 0;
					for(int unitIndex2 = k; unitIndex2 <= (int)expExponentialParameters.size(); unitIndex2++)
					{
						if(!ignore[unitIndex2-1]) sum += exponentialParameters[unitIndex2-1];
					}
					expNormalisingConstant(unitIndex-1, z-1) = exp(sum);
#ifndef NDEBUG
					assert(expNormalisingConstant(unitIndex-1, z-1) == expNormalisingConstant(unitIndex-1, z-1));
#endif
				}
				else
				{
					expNormalisingConstant(unitIndex-1, z-1) = expExponentialParameters[k-1] * expNormalisingConstant(unitIndex, z-2) + expNormalisingConstant(unitIndex, z-1);
#ifndef NDEBUG
					assert(expNormalisingConstant(unitIndex-1, z-1) == expNormalisingConstant(unitIndex-1, z-1));
#endif
				}
			}
			k--;
		}
	}
	void computeExponentialParameters(conditionalPoissonArgs& args)
	{
		std::vector<mpfr_class>& rescaledWeights = args.rescaledWeights;
		std::vector<mpfr_class>& weights = args.weights;
		int nUnits = (int)weights.size();

		args.exponentialParameters.resize(nUnits);
		args.expExponentialParameters.resize(nUnits);
		mpfr_class sumExponentialParameters = 0;
		int excluded = 0;
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i] && !args.zeroWeights[i])
			{
				args.expExponentialParameters[i] = rescaledWeights[i] / (1 - rescaledWeights[i]);
				args.exponentialParameters[i] = log(args.expExponentialParameters[i]);
				sumExponentialParameters += args.exponentialParameters[i];
			}
			else
			{
				excluded++;
				args.expExponentialParameters[i] = 0;
				args.exponentialParameters[i] = 0;
			}
		}
		//Rescale so the exponential parameters sum no zero
		mpfr_class toSubtract = sumExponentialParameters / (nUnits - excluded);
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i] && !args.zeroWeights[i])
			{
				args.exponentialParameters[i] -= toSubtract;
				args.expExponentialParameters[i] = exp(args.exponentialParameters[i]);
			}
		}

	}
}
