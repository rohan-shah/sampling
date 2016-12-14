#include "samplingBase.h"
namespace sampling
{
	void samplingBase(int n, std::vector<int>& indices, std::vector<double>& weights, std::vector<bool>& zeroWeights, std::vector<bool>& deterministicInclusion, int& nDeterministic, int& nZeroWeights)
	{
		indices.clear();
		int nUnits = (int)weights.size();

		//Work out which units have weights zero on one
		nZeroWeights = 0;
		nDeterministic = 0;
		zeroWeights.resize(nUnits);
		deterministicInclusion.resize(nUnits);
		std::fill(zeroWeights.begin(), zeroWeights.end(), false);
		std::fill(deterministicInclusion.begin(), deterministicInclusion.end(), false);
		for(int i = 0; i < nUnits; i++)
		{
			if(weights[i] > 1 || weights[i] < 0) throw std::runtime_error("Weights must be between 0 and 1"); 
			if(weights[i] == 0)
			{
				nZeroWeights++;
				zeroWeights[i] = true;
			}
			else if(weights[i] == 1)
			{
				nDeterministic++;
				deterministicInclusion[i] = true;
				indices.push_back(i);
			}
		}
		if(n > nUnits - nZeroWeights)
		{
			throw std::runtime_error("Cannot select more units than there are units with non-zero weights");
		}
		else if(nDeterministic > n)
		{
			throw std::runtime_error("Cannot select fewer units than there are units with weight 1");
		}
		if(n == nUnits - nZeroWeights)
		{
			nDeterministic = n;
			for(int i = 0; i < nUnits; i++)
			{
				if(!zeroWeights[i] && !deterministicInclusion[i])
				{
					deterministicInclusion[i] = true;
					indices.push_back(i);
				}
			}
			return;
		}
	}
	void sampfordBase(int n, std::vector<int>& indices, std::vector<double>& weights, std::vector<double>& rescaledWeights, std::vector<bool>& zeroWeights, std::vector<bool>& deterministicInclusion, int& nDeterministic, int& nZeroWeights)
	{
		indices.clear();
		int nUnits = (int)weights.size();
		rescaledWeights.resize(nUnits);

		//Work out which units have zero weights or deterministic inclusion
		nZeroWeights = 0;
		nDeterministic = 0;
		zeroWeights.resize(nUnits);
		deterministicInclusion.resize(nUnits);
		std::fill(zeroWeights.begin(), zeroWeights.end(), false);
		std::fill(deterministicInclusion.begin(), deterministicInclusion.end(), false);
		for(int i = 0; i < nUnits; i++)
		{
			if(weights[i] < 0) throw std::runtime_error("Weights must be between 0 and 1"); 
			if(weights[i] == 0)
			{
				nZeroWeights++;
				zeroWeights[i] = true;
			}
		}
		if(n > nUnits - nZeroWeights)
		{
			throw std::runtime_error("Cannot select more units than there are units with non-zero weights");
		}
		if(n == nUnits - nZeroWeights)
		{
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
			}
			return;
		}
		//Work out which units are going to be deterministically selected. 
		double cumulative;
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
			double maxAllowed =  cumulative / (n - indices.size());
			//Any weights that are too big are included with probability 1
			for(int i = 0; i < nUnits; i++)
			{
				if(weights[i] >= maxAllowed && !deterministicInclusion[i])
				{
					deterministicInclusion[i] = true;
					indices.push_back(i);
					hasDeterministic = true;
					nDeterministic++;
				}
		}
		} while(hasDeterministic);
		
		if(nDeterministic > n)
		{
			throw std::runtime_error("Cannot select fewer units than there are units with weight 1");
		}

		if(nDeterministic == (int)n)
		{
			if(nDeterministic + nZeroWeights != nUnits)
			{
				throw std::runtime_error("There were units with probability zero of selection but non-zero weights");
			}
			for(int i = 0; i < nUnits; i++)
			{
				if(deterministicInclusion[i])
				{
					rescaledWeights[i] = 1;
				}
				else
				{
					if(!zeroWeights[i]) throw std::runtime_error("Internal error");
					rescaledWeights[i] = 0;
				}
			}
			return;
		}

		//Rescale the weights so that they sum to n
		double factor = (n - nDeterministic) / cumulative;
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
	}
}

