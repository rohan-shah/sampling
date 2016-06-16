#include "samplingBase.h"
namespace sampling
{
	bool samplingBase(int n, std::vector<int>& indices, std::vector<mpfr_class>& weights, std::vector<mpfr_class>& rescaledWeights, std::vector<bool>& zeroWeights, std::vector<bool>& deterministicInclusion, std::vector<mpfr_class>& inclusionProbabilities)
	{
		indices.clear();
		int nUnits = (int)weights.size();
		rescaledWeights.resize(nUnits);
		inclusionProbabilities.resize(nUnits);

		//Work out which units have zero weights
		int nZeros = 0;
		zeroWeights.resize(nUnits);
		std::fill(zeroWeights.begin(), zeroWeights.end(), false);
		for(int i = 0; i < nUnits; i++)
		{
			if(weights[i] == 0)
			{
				nZeros++;
				zeroWeights[i] = true;
				inclusionProbabilities[i] = rescaledWeights[i] = 0;
			}
		}
		deterministicInclusion.resize(nUnits);
		if(n > nUnits - nZeros)
		{
			throw std::runtime_error("Input n was too big, or too many units had zero weights");
		}
		else if(n == nUnits - nZeros)
		{
			indices.reserve(nUnits);
			for(int i = 0; i < nUnits; i++)
			{
				if(!zeroWeights[i])
				{
					rescaledWeights[i] = inclusionProbabilities[i] = 1;
					deterministicInclusion[i] = true;
					indices.push_back(i);
				}
				else
				{
					rescaledWeights[i] = inclusionProbabilities[i] = 0;
					deterministicInclusion[i] = false;
				}
			}
			return true;
		}
		std::fill(deterministicInclusion.begin(), deterministicInclusion.end(), false);
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
				if(!deterministicInclusion[i]) cumulative += weights[i];
			}
			mpfr_class maxAllowed =  cumulative / mpfr_class(n - indices.size());
			//Any weights that are too big are included with probability 1
			for(int i = 0; i < nUnits; i++)
			{
				if(weights[i].convert_to<double>() >= maxAllowed.convert_to<double>() && !deterministicInclusion[i])
				{
					deterministicInclusion[i] = true;
					indices.push_back(i);
					hasDeterministic = true;
					inclusionProbabilities[i] = 1;
				}
		}
		} while(hasDeterministic);
		int deterministicIndices = (int)indices.size();

		if(deterministicIndices == (int)n)
		{
			if(deterministicIndices + nZeros != nUnits)
			{
				throw std::runtime_error("There were units with probability zero of selection but non-zero weights");
			}
			for(int i = 0; i < nUnits; i++)
			{
				if(deterministicInclusion[i])
				{
					inclusionProbabilities[i] = rescaledWeights[i] = 1;
				}
				else
				{
					if(!zeroWeights[i]) throw std::runtime_error("Internal error");
					inclusionProbabilities[i] = rescaledWeights[i] = 0;
				}
			}
			return true;
		}

		//Rescale the weights so that they sum to n
		mpfr_class factor = mpfr_class(n - deterministicIndices)/ cumulative;
		if(cumulative == 0)
		{
			throw std::runtime_error("Divide by zero encountered");
		}
		//And also work out the exponential parameters
		for(int i = 0; i < nUnits; i++)
		{
			if(deterministicInclusion[i])
			{
				rescaledWeights[i] = 1;
			}
			else if(zeroWeights[i])
			{
				rescaledWeights[i] = 0;
			}
			else rescaledWeights[i] = weights[i]*factor;
		}
		return false;
	}
}

