#include <boost/test/unit_test.hpp>
#include "conditionalPoissonSequential.h"
BOOST_AUTO_TEST_CASE(conditionalPoissonSequential1, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonSequentialArgs args;
	std::vector<int> indices;
	std::vector<sampling::mpfr_class> inclusionProbabilities, weights, rescaledWeights;
	args.indices = &indices;
	args.inclusionProbabilities = &inclusionProbabilities;
	args.weights = &weights;
	args.rescaledWeights = &rescaledWeights;

	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.push_back(1);
	weights.push_back(2);
	weights.push_back(3);
	args.n = 2;
	int table[2];
	table[0] = table[1] = 0;
	for(int i = 0; i < 10000; i++)
	{
		conditionalPoissonSequential(args, randomSource);
		std::sort(indices.begin(), indices.end());
		BOOST_TEST((int)indices.size() == 2);
		BOOST_TEST(indices[1] == 2);
		table[indices[0]]++;
	}
	BOOST_TEST(table[0]/10000.0 == 0.2, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[1]/10000.0 == 0.8, boost::test_tools::tolerance(0.02));
}
BOOST_AUTO_TEST_CASE(conditionalPoissonSequential2, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonSequentialArgs args;
	std::vector<int> indices;
	std::vector<sampling::mpfr_class> inclusionProbabilities, weights, rescaledWeights;
	args.indices = &indices;
	args.inclusionProbabilities = &inclusionProbabilities;
	args.weights = &weights;
	args.rescaledWeights = &rescaledWeights;

	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.push_back(1);
	weights.push_back(2);
	weights.push_back(2);
	args.n = 2;
	int table[3];
	table[0] = table[1] = table[2] = 0;
	for(int i = 0; i < 20000; i++)
	{
		conditionalPoissonSequential(args, randomSource);
		std::sort(indices.begin(), indices.end());
		BOOST_TEST((int)indices.size() == 2);
		if(indices[0] == 0 && indices[1] == 1) table[0]++;
		if(indices[0] == 0 && indices[1] == 2) table[1]++;
		if(indices[0] == 1 && indices[1] == 2) table[2]++;
	}
	BOOST_TEST(table[0]/20000.0 == 0.125, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[1]/20000.0 == 0.125, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[2]/20000.0 == 0.75, boost::test_tools::tolerance(0.02));
}
