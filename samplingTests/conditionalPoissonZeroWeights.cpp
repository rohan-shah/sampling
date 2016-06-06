#include <boost/test/unit_test.hpp>
#include "conditionalPoisson.h"
BOOST_AUTO_TEST_CASE(conditionalPoissonZeroWeights1, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonArgs args;
	std::vector<int> indices;
	std::vector<sampling::mpfr_class> inclusionProbabilities, weights, rescaledWeights;
	args.indices = &indices;
	args.inclusionProbabilities = &inclusionProbabilities;
	args.weights = &weights;
	args.rescaledWeights = &rescaledWeights;
	args.calculateInclusionProbabilities = true;

	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.push_back(1);
	weights.push_back(2);
	weights.push_back(0);
	weights.push_back(0);
	args.n = 1;
	for(int i = 0; i < 100; i++)
	{
		sampling::conditionalPoisson(args, randomSource);
		if(indices[0] == 0)
		{
			BOOST_TEST(inclusionProbabilities[0].convert_to<double>() == 1.0/5.0);
		}
		else
		{
			BOOST_TEST(inclusionProbabilities[1].convert_to<double>() == 4.0/5.0);
		}
		BOOST_TEST(rescaledWeights[0].convert_to<double>() == 1.0/3.0);
		BOOST_TEST(rescaledWeights[1].convert_to<double>() == 2.0/3.0);
		BOOST_TEST(!args.deterministicInclusion[0]);
		BOOST_TEST(!args.deterministicInclusion[1]);
		BOOST_TEST(!args.zeroWeights[0]);
		BOOST_TEST(!args.zeroWeights[1]);
		for(int j = 2; j < 4; j++)
		{
			BOOST_TEST(args.zeroWeights[j]);
			BOOST_TEST(rescaledWeights[j].convert_to<double>() == 0);
		}
	}
}
BOOST_AUTO_TEST_CASE(conditionalPoissonZeroWeights2, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonArgs args;
	std::vector<int> indices;
	std::vector<sampling::mpfr_class> inclusionProbabilities, weights, rescaledWeights;
	args.indices = &indices;
	args.inclusionProbabilities = &inclusionProbabilities;
	args.weights = &weights;
	args.rescaledWeights = &rescaledWeights;
	args.calculateInclusionProbabilities = true;

	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.push_back(0);
	weights.push_back(0);
	weights.push_back(0);
	weights.push_back(0);

	args.n = 1;
	BOOST_CHECK_THROW(sampling::conditionalPoisson(args, randomSource), std::runtime_error);
	args.n = 2;
	BOOST_CHECK_THROW(sampling::conditionalPoisson(args, randomSource), std::runtime_error);
}

