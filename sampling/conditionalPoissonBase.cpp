#include "conditionalPoissonBase.h"
namespace sampling
{
	conditionalPoissonArgs::conditionalPoissonArgs(conditionalPoissonArgs&& other)
		: indices(std::move(other.indices)), weights(std::move(other.weights)), n(other.n), expExponentialParameters(std::move(other.expExponentialParameters)), expNormalisingConstant(std::move(other.expNormalisingConstant)), deterministicInclusion(std::move(other.deterministicInclusion)), zeroWeights(std::move(other.zeroWeights))
	{}
	void calculateExpNormalisingConstants(conditionalPoissonArgs& args)
	{
		std::vector<mpfr_class>& weights = args.weights;
		std::vector<bool>& zeroWeights = args.zeroWeights;
		std::vector<bool>& deterministicInclusion = args.deterministicInclusion;
		int nUnits = (int)weights.size();
		//Now compute the inclusion probabilities
		calculateExpNormalisingConstants(args.expExponentialParameters, args.expNormalisingConstant, (int)args.n, nUnits, zeroWeights, deterministicInclusion);
	}
	void conditionalPoissonSecondOrderInclusionProbabilities(conditionalPoissonArgs& args, std::vector<mpfr_class>& inclusionProbabilities, boost::numeric::ublas::matrix<mpfr_class>& secondOrder)
	{
		std::size_t nUnits = inclusionProbabilities.size();
		secondOrder.resize(nUnits, nUnits);
		std::vector<int> identical;
		for(std::size_t i = 0; i < nUnits; i++)
		{
			mpfr_class remaining = args.n * inclusionProbabilities[i];
			identical.clear();
			for(std::size_t j = 0; j < nUnits; j++)
			{
				if(i == j)
				{
					secondOrder(i, j) = inclusionProbabilities[i];
					remaining -= secondOrder(i, j);
				}
				else if(i != j && args.expExponentialParameters[i].convert_to<double>() == args.expExponentialParameters[j].convert_to<double>())
				{
					identical.push_back((int)j);
				}
				else
				{
					secondOrder(i, j) = (args.expExponentialParameters[i]/(args.expExponentialParameters[i] - args.expExponentialParameters[j])) * inclusionProbabilities[j] + (args.expExponentialParameters[j]/(args.expExponentialParameters[j] - args.expExponentialParameters[i])) * inclusionProbabilities[i];
					remaining -= secondOrder(i, j);
				}
			}
			remaining /= identical.size();
			for(std::vector<int>::iterator k = identical.begin(); k != identical.end(); k++)
			{
				secondOrder(i, *k) = remaining;
			}
		}
	}
	void conditionalPoissonInclusionProbabilities(conditionalPoissonArgs& args, std::vector<mpfr_class>& inclusionProbabilities)
	{
		std::vector<mpfr_class>& weights = args.weights;
		std::vector<bool>& zeroWeights = args.zeroWeights;
		std::vector<bool>& deterministicInclusion = args.deterministicInclusion;
		int nUnits = (int)weights.size();
		//Now compute the inclusion probabilities
		inclusionProbabilities.resize(nUnits);
		int deterministicIndices = 0, nZeroWeights = 0;
		for(int i = 0; i < nUnits; i++)
		{
			if(deterministicInclusion[i])
			{
				deterministicIndices++;
				inclusionProbabilities[i] = 1;
			}
			else if(zeroWeights[i])
			{
				inclusionProbabilities[i] = 0;
				nZeroWeights++;
			}
		}
		if((int)args.n == deterministicIndices) return;
		else if((int)args.n < deterministicIndices)
		{
			throw std::runtime_error("Sample size is smaller than the number of deterministically included units");
		}
		else if(nUnits - nZeroWeights < (int)args.n)
		{
			throw std::runtime_error("Sample size is larger than the number of units with non-zero weights");
		}
		calculateExpNormalisingConstants(args.expExponentialParameters, args.expNormalisingConstant, (int)args.n, nUnits, zeroWeights, deterministicInclusion);
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
					inclusionProbabilities[unitCounter] -= args.expNormalisingConstant(0, j - 2) / boost::multiprecision::pow(args.expExponentialParameters[unitCounter], (j - 1));
				}
				else
				{
					inclusionProbabilities[unitCounter] += args.expNormalisingConstant(0, j - 2) / boost::multiprecision::pow(args.expExponentialParameters[unitCounter], (j - 1));
				}
			}
			inclusionProbabilities[unitCounter] *= boost::multiprecision::pow(args.expExponentialParameters[unitCounter], (args.n - deterministicIndices)) / expNormalisingConstant;
			if(inclusionProbabilities[unitCounter] <= 0)
			{
				std::stringstream ss;
				ss << "Inclusion probability had negative value " << inclusionProbabilities[unitCounter] << ", probably because of numerical instability";
				throw std::runtime_error(ss.str().c_str());
			}
		}
	}
	void calculateExpNormalisingConstants(std::vector<mpfr_class>& expExponentialParameters, boost::numeric::ublas::matrix<mpfr_class>& expNormalisingConstant, int n, int nUnits, std::vector<bool>& zeroWeights, std::vector<bool>& deterministicInclusion)
	{
		int nZeroWeights = 0, nDeterministic = 0;;
		for(int i = 0; i < nUnits; i++)
		{
			if(zeroWeights[i]) nZeroWeights++;
			if(deterministicInclusion[i]) nDeterministic++;
		}
		//We start by computing the normalising constants. First index is k, second is z. All indices in this loop are 1 indexed

		//Both constructing *and* destructing mpfr_class objects is expensive, so were need to make as few changes to the dimension here as possible. 
		if ((int)expNormalisingConstant.size1() < nUnits - nZeroWeights - nDeterministic || (int)expNormalisingConstant.size2() < n - nDeterministic)
		{
			expNormalisingConstant.resize(std::max(nUnits - nZeroWeights - nDeterministic, (int)expNormalisingConstant.size1()), std::max(n - nDeterministic, (int)expNormalisingConstant.size2()), false);
		}
		//This will skip over the *ignored* units (the ones that were deterministically included)
		int k = (int)expExponentialParameters.size();
		for(int unitIndex = nUnits - nDeterministic - nZeroWeights; unitIndex >= 1; unitIndex--)
		{
			while(zeroWeights[k-1] || deterministicInclusion[k-1]) k--;
			for(int z = 1; z <= std::min(n - nDeterministic, nUnits - nZeroWeights - nDeterministic - unitIndex+1); z++)
			{
				if(z == 1)
				{
					mpfr_class sum = 0;
					for(int unitIndex2 = k; unitIndex2 <= (int)expExponentialParameters.size(); unitIndex2++)
					{
						if(!zeroWeights[unitIndex2-1] && !deterministicInclusion[unitIndex2-1]) sum += expExponentialParameters[unitIndex2-1];
					}
					expNormalisingConstant(unitIndex-1, z-1) = sum;
#ifndef NDEBUG
					assert(expNormalisingConstant(unitIndex-1, z-1) == expNormalisingConstant(unitIndex-1, z-1));
#endif
				}
				else if(z == nUnits - nZeroWeights - nDeterministic - unitIndex + 1)
				{
					mpfr_class product = 1;
					for(int unitIndex2 = k; unitIndex2 <= (int)expExponentialParameters.size(); unitIndex2++)
					{
						if(!zeroWeights[unitIndex2-1] && !deterministicInclusion[unitIndex2-1]) product *= expExponentialParameters[unitIndex2-1];
					}
					expNormalisingConstant(unitIndex-1, z-1) = product;
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
		std::vector<mpfr_class>& weights = args.weights;
		int nUnits = (int)weights.size();

		args.expExponentialParameters.resize(nUnits);
		mpfr_class productExponentialParameters = 1;
		int excluded = 0;
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i] && !args.zeroWeights[i])
			{
				args.expExponentialParameters[i] = weights[i] / (1 - weights[i]);
				productExponentialParameters *= args.expExponentialParameters[i];
			}
			else
			{
				excluded++;
				args.expExponentialParameters[i] = 0;
			}
		}
		//Rescale so the exponential parameters sum no zero
		mpfr_class expToDivide = boost::multiprecision::pow(productExponentialParameters, 1.0 / (nUnits - excluded));
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i] && !args.zeroWeights[i])
			{
				args.expExponentialParameters[i] /= expToDivide;
			}
		}
	}
}
