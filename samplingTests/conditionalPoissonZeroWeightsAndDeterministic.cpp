#include <boost/test/unit_test.hpp>
#include "conditionalPoissonRejective.h"
BOOST_AUTO_TEST_CASE(conditionalPoissonRejectiveZeroWeightsAndDeterministic1, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonRejectiveArgs args(true);
	std::vector<int>& indices = args.indices;
	std::vector<sampling::mpfr_class>& inclusionProbabilities = args.inclusionProbabilities;
	std::vector<sampling::mpfr_class>& weights = args.weights;
	std::vector<bool>& zeroWeights = args.zeroWeights;
	std::vector<bool>& deterministicInclusion = args.deterministicInclusion;

	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.clear();
	weights.push_back(1.0/3.0);
	weights.push_back(2.0/3.0);
	weights.push_back(0);
	weights.push_back(1);
	args.n = 2;
	for(int i = 0; i < 100; i++)
	{
		sampling::conditionalPoissonRejective(args, randomSource);
		std::sort(indices.begin(), indices.end());
		BOOST_TEST(indices.size() == (std::size_t)2);
		BOOST_TEST(inclusionProbabilities.size() == (std::size_t)4);
		BOOST_TEST(weights.size() == (std::size_t)4);
		BOOST_TEST(inclusionProbabilities[2].convert_to<double>() == 0.0);
		BOOST_TEST(inclusionProbabilities[3].convert_to<double>() == 1.0);
		BOOST_TEST(zeroWeights[2]);
		BOOST_TEST(deterministicInclusion[3]);
		BOOST_TEST(indices[1] == 3);
		BOOST_TEST((indices[0] == 0 || indices[0] == 1));
		if(indices[0] == 0)
		{
			BOOST_TEST(inclusionProbabilities[0].convert_to<double>() == 1.0/5.0);
		}
		else
		{
			BOOST_TEST(inclusionProbabilities[1].convert_to<double>() == 4.0/5.0);
		}
	}	
}
BOOST_AUTO_TEST_CASE(conditionalPoissonRejectiveZeroWeightsAndDeterministic2, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonRejectiveArgs args(true);
	std::vector<int>& indices = args.indices;
	std::vector<sampling::mpfr_class>& inclusionProbabilities = args.inclusionProbabilities;
	std::vector<sampling::mpfr_class>& weights = args.weights;
	
	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.push_back(0);
	weights.push_back(0);
	weights.push_back(1.0);
	weights.push_back(1.0);

	args.n = 2;
	for(int i = 0; i < 100; i++)
	{
		sampling::conditionalPoissonRejective(args, randomSource);
		std::sort(indices.begin(), indices.end());
		BOOST_TEST(indices.size() == (std::size_t)2);
		BOOST_TEST(inclusionProbabilities.size() == (std::size_t)4);
		BOOST_TEST(weights.size() == (std::size_t)4);
		BOOST_TEST(inclusionProbabilities[0].convert_to<double>() == 0.0);
		BOOST_TEST(inclusionProbabilities[1].convert_to<double>() == 0.0);
		BOOST_TEST(inclusionProbabilities[2].convert_to<double>() == 1.0);
		BOOST_TEST(inclusionProbabilities[3].convert_to<double>() == 1.0);
		BOOST_TEST(indices[0] == 2);
		BOOST_TEST(indices[1] == 3);
	}
}

