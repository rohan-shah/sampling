#include <boost/test/unit_test.hpp>
#include "conditionalPoissonDrafting.h"
BOOST_AUTO_TEST_CASE(conditionalPoissonDrafting1, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonDraftingArgs args;
	std::vector<int>& indices = args.indices;
	std::vector<sampling::mpfr_class>& weights = args.weights;

	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.push_back(1);
	weights.push_back(2);
	weights.push_back(3);
	args.n = 2;
	int table[2];
	table[0] = table[1] = 0;
	for(int i = 0; i < 5000; i++)
	{
		conditionalPoissonDrafting(args, randomSource);
		std::sort(indices.begin(), indices.end());
		BOOST_TEST((int)indices.size() == 2);
		BOOST_TEST(indices[1] == 2);
		table[indices[0]]++;
	}
	BOOST_TEST(table[0]/5000.0 == 0.2, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[1]/5000.0 == 0.8, boost::test_tools::tolerance(0.02));
}
BOOST_AUTO_TEST_CASE(conditionalPoissonDrafting2, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonDraftingArgs args;
	std::vector<int>& indices = args.indices;
	std::vector<sampling::mpfr_class>& weights = args.weights;

	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.push_back(1);
	weights.push_back(2);
	weights.push_back(2);
	args.n = 2;
	int table[3];
	table[0] = table[1] = table[2] = 0;
	for(int i = 0; i < 5000; i++)
	{
		conditionalPoissonDrafting(args, randomSource);
		std::sort(indices.begin(), indices.end());
		BOOST_TEST((int)indices.size() == 2);
		if(indices[0] == 0 && indices[1] == 1) table[0]++;
		if(indices[0] == 0 && indices[1] == 2) table[1]++;
		if(indices[0] == 1 && indices[1] == 2) table[2]++;
	}
	BOOST_TEST(table[0]/5000.0 == 0.125, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[1]/5000.0 == 0.125, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[2]/5000.0 == 0.75, boost::test_tools::tolerance(0.02));
}
