#include "systematicSampling.h"
#include <boost/random/uniform_real_distribution.hpp>
namespace sampling
{
	void systematicSampling(const std::vector<mpfr_class>& weights, mpfr_class interval, std::vector<int>& indices, boost::mt19937& randomSource)
	{
		indices.clear();
		double position = boost::random::uniform_real_distribution<>(0, interval.convert_to<double>())(randomSource);
		for(int unit = 0; unit < (int)weights.size(); unit++)
		{
			position -= weights[unit].convert_to<double>();
			if(position < 0)
			{
				indices.push_back(unit);
				position += interval.convert_to<double>();
			}
		}
	}
}
namespace samplingDouble
{
	void systematicSamplingDouble(const std::vector<double>& weights, double interval, std::vector<int>& indices, boost::mt19937& randomSource)
	{
		indices.clear();
		double position = boost::random::uniform_real_distribution<>(0, interval)(randomSource);
		for(int unit = 0; unit < (int)weights.size(); unit++)
		{
			position -= weights[unit];
			if(position < 0)
			{
				indices.push_back(unit);
				position += interval;
			}
		}
	}
}
