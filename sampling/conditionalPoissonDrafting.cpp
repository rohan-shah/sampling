#include "conditionalPoissonDrafting.h"
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "samplingBase.h"
namespace sampling
{
	using boost::multiprecision::log;
	using std::log;
	using boost::multiprecision::exp;
	using std::exp;
	void conditionalPoissonDrafting(conditionalPoissonDraftingArgs& args, boost::mt19937& randomSource)
	{
		std::vector<int>& indices = *args.indices;
		std::vector<mpfr_class>& weights = *args.weights;
		std::vector<mpfr_class>& inclusionProbabilities = *args.inclusionProbabilities;
		std::vector<mpfr_class>& rescaledWeights = *args.rescaledWeights;
		int nUnits = (int)weights.size();

		if(samplingBase(args.n, indices, weights, rescaledWeights, args.zeroWeights, args.deterministicInclusion, inclusionProbabilities)) return;
		int deterministicIndices = indices.size();

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
		conditionalPoissonInclusionProbabilities(args);
		args.inclusionProbabilities2.clear();
		args.inclusionProbabilities2.insert(args.inclusionProbabilities2.begin(), inclusionProbabilities.begin(), inclusionProbabilities.end());
		//drafting procedure
		std::vector<int> remaining;
		remaining.reserve(nUnits);
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i] && ! args.zeroWeights[i]) remaining.push_back(i);
		}
		std::vector<double> accumulated;
		accumulated.reserve(nUnits);
		boost::random::uniform_real_distribution<double> uniformDest;
		for(int i = 0; i < args.n - deterministicIndices; i++)
		{
			accumulated.clear();
			mpfr_class sum = 0;
			//Accumulate inclusion probabilities
			for(int j = 0; j < remaining.size(); j++)
			{
				sum += args.inclusionProbabilities2[remaining[j]];
				accumulated.push_back(sum.convert_to<double>());
			}
			//Sample according to inclusion probabilities
			double uniformSample = uniformDest(randomSource);
			int selectedUnit = -1, selectedUnitInRemaining = -1;
			for(int j = 0; j < accumulated.size(); j++)
			{
				if(uniformSample <= accumulated[j])
				{
					selectedUnit = remaining[j];
					selectedUnitInRemaining = j;
				}
			}
			if(selectedUnit == -1) throw std::runtime_error("Internal error");
			//Update the inclusion probabilities to reflect the sampled unit. See Rare Event Simulation Using Monte Carlo Methods, edited by Gerardo Rubino and Bruno Tuffin, p183. 
			args.inclusionProbabilities3.resize(nUnits);
			std::swap(remaining.back(), remaining[selectedUnitInRemaining]);
			remaining.pop_back();
			if(args.n - deterministicIndices - i - 1 == 0) throw std::runtime_error("Internal error");
			for(int j = 0; j < remaining.size(); j++)
			{
				args.inclusionProbabilities3[remaining[j]] = (args.expExponentialParameters[selectedUnit] * args.inclusionProbabilities2[remaining[j]] - args.expExponentialParameters[remaining[j]] * args.inclusionProbabilities2[selectedUnit]) / ((args.n - deterministicIndices - i - 1) * (args.expExponentialParameters[selectedUnit] - args.expExponentialParameters[remaining[j]]) * args.inclusionProbabilities2[selectedUnit]);
			}
			args.inclusionProbabilities3.swap(args.inclusionProbabilities2);
		}
	}
}
