#include "fearnheadGetKappa.h"
namespace samplingDouble
{
	void fearnheadGetKappa(std::vector<double>& sortedWeights, boost::mt19937& randomSource, int N, int& A, double& B)
	{
#ifndef NDEBUG
		if(!std::is_sorted(sortedWeights.begin(), sortedWeights.end()))
		{
			throw std::runtime_error("Internal error");
		}
#endif
		std::size_t nUnits = sortedWeights.size();
		int lowerBound = 0, upperBound = (int)nUnits-1;
		double lowerBoundSum = 0;
		const double tolerance = 1e-8;
		while(true)
		{
			//Split the existing partition between currentIndex and currentIndex2
			int currentIndex = (lowerBound + upperBound - 1) / 2;
			while(currentIndex <= upperBound - 1 && sortedWeights[currentIndex] == sortedWeights[currentIndex+1]) currentIndex++;
			if(currentIndex == (int)nUnits - 1)
			{
				double partitionWeight = sortedWeights[currentIndex];
				A = (int)nUnits - lowerBound;
				B = lowerBoundSum;
				int i = lowerBound;
				for(; i < (int)nUnits && sortedWeights[i] < partitionWeight; i++)
				{
					A--;
					B += sortedWeights[i];
				}
				if(B / partitionWeight + A > N + tolerance)
				{
					A = 0;
					for(; i < (int)nUnits; i++)
					{
						B += sortedWeights[i];
					}
					return;
				}
				else upperBound = (lowerBound + upperBound - 1) / 2;
			}
			else
			{
				double partitionWeight1 = sortedWeights[currentIndex], partitionWeight2 = sortedWeights[currentIndex+1];
				double B1 = lowerBoundSum;
				int A1 = (int)nUnits - lowerBound;
				int i = lowerBound;
				for(; sortedWeights[i] < partitionWeight1 && i < (int)nUnits; i++)
				{
					B1 += sortedWeights[i];
					A1--;
				}
				double B2 = B1;
				int A2 = A1;
				for(; sortedWeights[i] < partitionWeight2 && i < (int)nUnits; i++)
				{
					B2 += sortedWeights[i];
					A2--;
				}
				//We've found the minimum
				if(B1 / partitionWeight1 + A1 > N + tolerance && B2 / partitionWeight2 + A2 <= N + tolerance)
				{
					B = B2;
					A = A2;
					return;
				}
				//Take the smaller half
				else if(B1 / partitionWeight1 + A1 <= N + tolerance)
				{
					upperBound = (lowerBound + upperBound - 1) / 2;
				}
				//Take the bigger half
				else
				{
					lowerBound = currentIndex + 1;
					lowerBoundSum = B2;
				}
			}
			if(sortedWeights[lowerBound] == sortedWeights[upperBound])
			{
				double partitionWeight = sortedWeights[lowerBound];
				B = lowerBoundSum;
				A = (int)nUnits - lowerBound;
				int i = lowerBound;
				for(; sortedWeights[i] < partitionWeight; i++)
				{
					B += sortedWeights[i];
					A--;
				}
				if(B / partitionWeight + A <= N + tolerance)
				{
					return;
				}
				else
				{
					for(; i <= upperBound; i++)
					{
						B += sortedWeights[i];
						A--;
					}
					return;
				}
			}
		}
	}
}
