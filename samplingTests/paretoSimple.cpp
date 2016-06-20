#include <boost/test/unit_test.hpp>
#include "pareto.h"
BOOST_AUTO_TEST_CASE(paretoSimple1, * boost::unit_test::tolerance(0.00001))
{
	sampling::paretoSamplingArgs args;
	std::vector<int>& indices = args.indices;
	std::vector<sampling::mpfr_class>& weights = args.weights;

	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.push_back(1.0/3.0);
	weights.push_back(2.0/3.0);
	args.n = 1;
	for(int i = 0; i < 100; i++)
	{
		sampling::pareto(args, randomSource);
		BOOST_TEST(!args.deterministicInclusion[0]);
		BOOST_TEST(!args.deterministicInclusion[1]);
		BOOST_TEST(!args.zeroWeights[0]);
		BOOST_TEST(!args.zeroWeights[1]);
		BOOST_TEST(args.zeroWeights.size() == (std::size_t)2);
		BOOST_TEST(args.deterministicInclusion.size() == (std::size_t)2);
		BOOST_TEST(indices.size() == (std::size_t)1);
		BOOST_TEST(weights.size() == (std::size_t)2);
	}
}
BOOST_AUTO_TEST_CASE(paretoSimple2, * boost::unit_test::tolerance(0.00001))
{
	sampling::paretoSamplingArgs args;
	std::vector<int>& indices = args.indices;
	std::vector<sampling::mpfr_class>& weights = args.weights;

	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.push_back(3.0/6.0);
	weights.push_back(4.0/6.0);
	weights.push_back(5.0/6.0);
	args.n = 2;
	for(int i = 0; i < 100; i++)
	{
		sampling::pareto(args, randomSource);
		std::sort(indices.begin(), indices.end());
		BOOST_TEST(!args.deterministicInclusion[0]);
		BOOST_TEST(!args.deterministicInclusion[1]);
		BOOST_TEST(!args.deterministicInclusion[2]);
		BOOST_TEST(!args.zeroWeights[0]);
		BOOST_TEST(!args.zeroWeights[1]);
		BOOST_TEST(!args.zeroWeights[2]);
		BOOST_TEST(args.zeroWeights.size() == (std::size_t)3);
		BOOST_TEST(args.deterministicInclusion.size() == (std::size_t)3);
		BOOST_TEST(indices.size() == (std::size_t)2);
		BOOST_TEST(weights.size() == (std::size_t)3);
	}
}
