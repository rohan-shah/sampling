#include "conditionalPoisson.h"
#include <boost/random/bernoulli_distribution.hpp>
namespace sampling
{
	using boost::multiprecision::log;
	using std::log;
	using boost::multiprecision::exp;
	using std::exp;
	void conditionalPoisson(conditionalPoissonArgs& args, boost::mt19937& randomSource)
	{
		std::vector<int>& indices = *args.indices;
		std::vector<mpfr_class>& weights = *args.weights;
		std::vector<mpfr_class>& inclusionProbabilities = *args.inclusionProbabilities;
		std::vector<mpfr_class>& rescaledWeights = *args.rescaledWeights;

		int nUnits = (int)weights.size();

		indices.clear();
		rescaledWeights.resize(nUnits);
		inclusionProbabilities.resize(nUnits);

		//Work out which units have zero weights
		int nZeros = 0;
		args.zeroWeights.resize(nUnits);
		std::fill(args.zeroWeights.begin(), args.zeroWeights.end(), false);
		for(int i = 0; i < nUnits; i++)
		{
			if(weights[i] == 0)
			{
				nZeros++;
				args.zeroWeights[i] = true;
				inclusionProbabilities[i] = rescaledWeights[i] = 0;
			}
		}
		args.deterministicInclusion.resize(nUnits);
		if((int)args.n > nUnits - nZeros)
		{
			throw std::runtime_error("Input n was too big, or too many units had zero weights");
		}
		else if((int)args.n == nUnits - nZeros)
		{
			indices.reserve(nUnits);
			for(int i = 0; i < nUnits; i++)
			{
				if(!args.zeroWeights[i])
				{
					rescaledWeights[i] = inclusionProbabilities[i] = 1;
					args.deterministicInclusion[i] = true;
					indices.push_back(i);
				}
				else
				{
					rescaledWeights[i] = inclusionProbabilities[i] = 0;
					args.deterministicInclusion[i] = false;
				}
			}
			return;
		}
		std::fill(args.deterministicInclusion.begin(), args.deterministicInclusion.end(), false);
		//Work out which units are going to be deterministically selected. 
		mpfr_class cumulative;
		bool hasDeterministic = false;
		do
		{
			hasDeterministic = false;
			//Work out sum of weights
			cumulative = 0;
			for(int i = 0; i < nUnits; i++)
			{
				if(!args.deterministicInclusion[i]) cumulative += weights[i];
			}
			mpfr_class maxAllowed =  cumulative / mpfr_class(args.n - indices.size());
			//Any weights that are too big are included with probability 1
			for(int i = 0; i < nUnits; i++)
			{
				if(weights[i] >= maxAllowed && !args.deterministicInclusion[i])
				{
					args.deterministicInclusion[i] = true;
					indices.push_back(i);
					hasDeterministic = true;
					inclusionProbabilities[i] = 1;
				}
		}
		} while(hasDeterministic);
		int deterministicIndices = (int)indices.size();

		if(deterministicIndices == (int)args.n)
		{
			for(int i = 0; i < nUnits; i++)
			{
				if(args.deterministicInclusion[i])
				{
					inclusionProbabilities[i] = rescaledWeights[i] = 1;
				}
				else
				{
					if(!args.zeroWeights[i]) throw std::runtime_error("Internal error");
					inclusionProbabilities[i] = rescaledWeights[i] = 0;
				}
			}
			return;
		}

		//Rescale the weights so that they sum to n
		mpfr_class factor = mpfr_class(args.n - deterministicIndices)/ cumulative;
		if(cumulative == 0)
		{
			throw std::runtime_error("Divide by zero encountered");
		}
		//And also work out the exponential parameters
		for(int i = 0; i < nUnits; i++)
		{
			if(args.deterministicInclusion[i])
			{
				rescaledWeights[i] = 1;
			}
			else if(args.zeroWeights[i])
			{
				rescaledWeights[i] = 0;
			}
			else rescaledWeights[i] = weights[i]*factor;
		}

		args.exponentialParameters.resize(nUnits);
		args.expExponentialParameters.resize(nUnits);
		mpfr_class sumExponentialParameters = 0;
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
				args.expExponentialParameters[i] = 0;
				args.exponentialParameters[i] = 0;
			}
		}
		//Rescale so the exponential parameters sum no zero
		mpfr_class toSubtract = sumExponentialParameters / nUnits;
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i])
			{
				args.exponentialParameters[i] -= toSubtract;
				args.expExponentialParameters[i] = exp(args.exponentialParameters[i]);
			}
		}
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
	void conditionalPoissonInclusionProbabilities(conditionalPoissonArgs& args)
	{
		std::vector<mpfr_class>& inclusionProbabilities = *args.inclusionProbabilities;
		std::vector<mpfr_class>& rescaledWeights = *args.rescaledWeights;
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
				}
				else if(z == nUnits - unitIndex + 1)
				{
					mpfr_class sum = 0;
					for(int unitIndex2 = k; unitIndex2 <= (int)expExponentialParameters.size(); unitIndex2++)
					{
						if(!ignore[unitIndex2-1]) sum += exponentialParameters[unitIndex2-1];
					}
					expNormalisingConstant(unitIndex-1, z-1) = exp(sum);
				}
				else
				{
					expNormalisingConstant(unitIndex-1, z-1) = expExponentialParameters[k-1] * expNormalisingConstant(unitIndex, z-2) + expNormalisingConstant(unitIndex, z-1);
				}
			}
			k--;
		}
	}
}
