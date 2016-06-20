#include <boost/test/unit_test.hpp>
#include "conditionalPoissonRejective.h"
BOOST_AUTO_TEST_CASE(conditionalPoissonRejectiveZeroWeights1, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonRejectiveArgs args(true);
	std::vector<int>& indices = args.indices;
	std::vector<sampling::mpfr_class>& inclusionProbabilities = args.inclusionProbabilities;
	std::vector<sampling::mpfr_class>& weights = args.weights;

	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.push_back(1.0/3.0);
	weights.push_back(2.0/3.0);
	weights.push_back(0);
	weights.push_back(0);
	args.n = 1;
	for(int i = 0; i < 100; i++)
	{
		sampling::conditionalPoissonRejective(args, randomSource);
		if(indices[0] == 0)
		{
			BOOST_TEST(inclusionProbabilities[0].convert_to<double>() == 1.0/5.0);
		}
		else
		{
			BOOST_TEST(inclusionProbabilities[1].convert_to<double>() == 4.0/5.0);
		}
		BOOST_TEST(!args.deterministicInclusion[0]);
		BOOST_TEST(!args.deterministicInclusion[1]);
		BOOST_TEST(!args.zeroWeights[0]);
		BOOST_TEST(!args.zeroWeights[1]);
		for(int j = 2; j < 4; j++)
		{
			BOOST_TEST(args.zeroWeights[j]);
			BOOST_TEST(inclusionProbabilities[j].convert_to<double>() == 0);
		}
		BOOST_TEST(args.zeroWeights.size() == (std::size_t)4);
		BOOST_TEST(args.deterministicInclusion.size() == (std::size_t)4);
		BOOST_TEST(indices.size() == (std::size_t)1);
		BOOST_TEST(weights.size() == (std::size_t)4);
		BOOST_TEST(inclusionProbabilities.size() == (std::size_t)4);
	}
}
BOOST_AUTO_TEST_CASE(conditionalPoissonRejectiveZeroWeights2, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonRejectiveArgs args(true);
	std::vector<sampling::mpfr_class>& weights = args.weights;

	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.push_back(0);
	weights.push_back(0);
	weights.push_back(0);
	weights.push_back(0);

	args.n = 1;
	BOOST_CHECK_THROW(sampling::conditionalPoissonRejective(args, randomSource), std::runtime_error);
	args.n = 2;
	BOOST_CHECK_THROW(sampling::conditionalPoissonRejective(args, randomSource), std::runtime_error);
}

