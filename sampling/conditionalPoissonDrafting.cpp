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
		computeExponentialParameters(args);

		conditionalPoissonInclusionProbabilities(args, inclusionProbabilities);
		args.inclusionProbabilities2.clear();
		args.inclusionProbabilities2.insert(args.inclusionProbabilities2.begin(), inclusionProbabilities.begin(), inclusionProbabilities.end());
		//They have to be rescaled first
		for(int i = 0; i < (int)args.inclusionProbabilities2.size(); i++) args.inclusionProbabilities2[i] /= args.n - deterministicIndices;
		//drafting procedure
		std::vector<int> remaining;
		remaining.reserve(nUnits);
		for(int i = 0; i < nUnits; i++)
		{
			if(!args.deterministicInclusion[i] && ! args.zeroWeights[i]) remaining.push_back(i);
		}
		std::vector<double> accumulated;
		accumulated.reserve(nUnits);
		std::vector<int> equalProbabilityUnits;
		boost::random::uniform_real_distribution<double> uniformDest;
		for(int i = 0; i < (int)args.n - deterministicIndices; i++)
		{
			accumulated.clear();
			mpfr_class sum = 0;
			//Accumulate inclusion probabilities
			for(int j = 0; j < (int)remaining.size(); j++)
			{
				sum += args.inclusionProbabilities2[remaining[j]];
				accumulated.push_back(sum.convert_to<double>());
			}
			//Sample according to inclusion probabilities
			double uniformSample = uniformDest(randomSource);
			int selectedUnit = -1, selectedUnitInRemaining = -1;
			for(int j = 0; j < (int)accumulated.size(); j++)
			{
				if(uniformSample <= accumulated[j])
				{
					selectedUnit = remaining[j];
					selectedUnitInRemaining = j;
					break;
				}
			}
			if(selectedUnit == -1) throw std::runtime_error("Internal error");
			indices.push_back(selectedUnit);
			//Update the inclusion probabilities to reflect the sampled unit. See Rare Event Simulation Using Monte Carlo Methods, edited by Gerardo Rubino and Bruno Tuffin, p183. 
			args.inclusionProbabilities3.resize(nUnits);
			if(args.n - deterministicIndices - i - 1 == 0) 
			{
				break;
			}
			//Things get very messy in the case that there are multiple units with the same weight. See Sampling Algorithms by Yves Tille, pp 102. 
			mpfr_class sumNotEqual = 1;
			equalProbabilityUnits.clear();
			for(int j = 0; j < (int)remaining.size(); j++)
			{
				if(j != selectedUnitInRemaining)
				{
					if(args.expExponentialParameters[selectedUnit] == args.expExponentialParameters[remaining[j]])
					{
						equalProbabilityUnits.push_back(remaining[j]);
					}
					else
					{
						args.inclusionProbabilities3[remaining[j]] = (args.expExponentialParameters[selectedUnit] * args.inclusionProbabilities2[remaining[j]] - args.expExponentialParameters[remaining[j]] * args.inclusionProbabilities2[selectedUnit]) / ((args.n - deterministicIndices - i - 1) * (args.expExponentialParameters[selectedUnit] - args.expExponentialParameters[remaining[j]]) * args.inclusionProbabilities2[selectedUnit]);
						sumNotEqual -= args.inclusionProbabilities3[remaining[j]];
					}
				}
			}
			for(std::vector<int>::iterator j = equalProbabilityUnits.begin(); j != equalProbabilityUnits.end(); j++)
			{
				args.inclusionProbabilities3[*j] = sumNotEqual / equalProbabilityUnits.size();
			}
			std::swap(remaining.back(), remaining[selectedUnitInRemaining]);
			remaining.pop_back();
			args.inclusionProbabilities3.swap(args.inclusionProbabilities2);
		}
	}
}
