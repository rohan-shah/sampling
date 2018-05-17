#include "fearnheadSampling.h"
#include "fearnheadGetKappa.h"
#include <boost/random/uniform_real_distribution.hpp>
namespace samplingDouble
{
	void fearnheadSampling(fearnheadSamplingArgs& args, boost::mt19937& randomSource)
	{
		args.sortedWeights = args.weights;
		std::sort(args.sortedWeights.begin(), args.sortedWeights.end());
		int A;
		double B;
		fearnheadGetKappa(args.sortedWeights, randomSource, args.n, A, B);
		args.c = ((int)args.n - A) / B;
#ifndef NDEBUG
		double checkSum = 0;
		for(std::vector<double>::iterator i = args.sortedWeights.begin(); i != args.sortedWeights.end(); i++)
		{
			checkSum += std::min(*i * args.c, 1.0);
		}
		if(fabs(checkSum - args.n) > 1e-6)
		{
			throw std::runtime_error("Internal error");
		}
#endif
		double cInverse = 1/args.c;

		args.indices.clear();
		args.deterministicInclusion.resize(args.weights.size());
		std::fill(args.deterministicInclusion.begin(), args.deterministicInclusion.end(), false);

		double sumWeights = 0;
		for(int unit = 0; unit < (int)args.weights.size(); unit++)
		{
			double currentWeight = args.weights[unit];
			if(currentWeight >= cInverse)
			{
				args.indices.push_back(unit);
				args.deterministicInclusion[unit] = true;
			}
			else sumWeights += currentWeight;
		}

		double interval = sumWeights / (args.n - args.indices.size());
		if (interval <= 0) throw std::runtime_error("Internal error in fearnheadSampling");
		double position = boost::random::uniform_real_distribution<>(0, interval)(randomSource);
		for(int unit = 0; unit < (int)args.weights.size(); unit++)
		{
			if(args.weights[unit] < cInverse)
			{
				position -= args.weights[unit];
				if(position < 0)
				{
					args.indices.push_back(unit);
					position += interval;
				}
			}
		}
	}
}
namespace sampling
{
	void fearnheadSampling(fearnheadSamplingArgs& args, boost::mt19937& randomSource)
	{
		args.sortedWeights = args.weights;
		std::sort(args.sortedWeights.begin(), args.sortedWeights.end());
		int A;
		mpfr_class B;
		fearnheadGetKappa(args.sortedWeights, randomSource, args.n, A, B);
		args.c = ((int)args.n - A) / B;
#ifndef NDEBUG
		double checkSum = 0;
		for(std::vector<mpfr_class>::iterator i = args.sortedWeights.begin(); i != args.sortedWeights.end(); i++)
		{
			checkSum += std::min(mpfr_class(*i * args.c).convert_to<double>(), 1.0);
		}
		if(fabs(checkSum - args.n) > 1e-6)
		{
			throw std::runtime_error("Internal error");
		}
#endif
		mpfr_class cInverse = 1.0/args.c;

		args.indices.clear();
		args.deterministicInclusion.resize(args.weights.size());
		std::fill(args.deterministicInclusion.begin(), args.deterministicInclusion.end(), false);

		mpfr_class sumWeights = 0;
		for(int unit = 0; unit < (int)args.weights.size(); unit++)
		{
			mpfr_class currentWeight = args.weights[unit];
			if(currentWeight >= cInverse)
			{
				args.indices.push_back(unit);
				args.deterministicInclusion[unit] = true;
			}
			else sumWeights += currentWeight;
		}

		double interval = mpfr_class(sumWeights / (args.n - args.indices.size())).convert_to<double>();
		mpfr_class position = boost::random::uniform_real_distribution<>(0, interval)(randomSource);
		for(int unit = 0; unit < (int)args.weights.size(); unit++)
		{
			if(args.weights[unit] < cInverse)
			{
				position -= args.weights[unit];
				if(position < 0)
				{
					args.indices.push_back(unit);
					position += interval;
				}
			}
		}
	}
}
