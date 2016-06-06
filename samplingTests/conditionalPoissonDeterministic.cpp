#include <boost/test/unit_test.hpp>
#include "conditionalPoisson.h"
BOOST_AUTO_TEST_CASE(conditionalPossonDeterministic1, * boost::unit_test::tolerance(0.00001))
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
	weights.push_back(3);
	args.n = 2;
	for(int i = 0; i < 100; i++)
	{
		sampling::conditionalPoisson(args, randomSource);
		std::sort(indices.begin(), indices.end());
		BOOST_TEST((std::find(indices.begin(), indices.end(), 2) != indices.end()));
		BOOST_TEST(!args.deterministicInclusion[0]);
		BOOST_TEST(!args.deterministicInclusion[1]);
		BOOST_TEST(args.deterministicInclusion[2]);
		BOOST_TEST(inclusionProbabilities[2].convert_to<double>() == 1);
		BOOST_TEST(rescaledWeights[2].convert_to<double>() == 1);
		BOOST_TEST(!args.zeroWeights[0]);
		BOOST_TEST(!args.zeroWeights[1]);
		BOOST_TEST(!args.zeroWeights[2]);
	}
}
BOOST_AUTO_TEST_CASE(conditionalPossonDeterministic2, * boost::unit_test::tolerance(0.00001))
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
	weights.push_back(3);
	args.n = 3;
	for(int i = 0; i < 10; i++)
	{
		sampling::conditionalPoisson(args, randomSource);
		BOOST_TEST(args.deterministicInclusion[0]);
		BOOST_TEST(args.deterministicInclusion[1]);
		BOOST_TEST(args.deterministicInclusion[2]);
		BOOST_TEST(inclusionProbabilities[0].convert_to<double>() == 1);
		BOOST_TEST(inclusionProbabilities[1].convert_to<double>() == 1);
		BOOST_TEST(inclusionProbabilities[2].convert_to<double>() == 1);
		BOOST_TEST(rescaledWeights[0].convert_to<double>() == 1);
		BOOST_TEST(rescaledWeights[1].convert_to<double>() == 1);
		BOOST_TEST(rescaledWeights[2].convert_to<double>() == 1);
		BOOST_TEST(!args.zeroWeights[0]);
		BOOST_TEST(!args.zeroWeights[1]);
		BOOST_TEST(!args.zeroWeights[2]);
	}
}
BOOST_AUTO_TEST_CASE(conditionalPossonDeterministic3, * boost::unit_test::tolerance(0.00001))
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
	weights.push_back(3);
	for(int i = 0; i < 7; i++) weights.push_back(100);
	args.n = 8;
	for(int i = 0; i < 10; i++)
	{
		sampling::conditionalPoisson(args, randomSource);
		for(int j = 3; j < 10; j++) 
		{
			BOOST_TEST(args.deterministicInclusion[j]);
			BOOST_TEST(inclusionProbabilities[j].convert_to<double>() == 1);
			BOOST_TEST(rescaledWeights[j].convert_to<double>() == 1);
		}
		for(int j = 0; j < 10; j++)
		{
			BOOST_TEST(!args.zeroWeights[j]);
		}
		BOOST_TEST(rescaledWeights[0].convert_to<double>() == 1.0/6.0);
		BOOST_TEST(rescaledWeights[1].convert_to<double>() == 2.0/6.0);
		BOOST_TEST(rescaledWeights[2].convert_to<double>() == 3.0/6.0);
	}
}
