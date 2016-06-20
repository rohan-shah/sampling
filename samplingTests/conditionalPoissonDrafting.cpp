#include <boost/test/unit_test.hpp>
#include "conditionalPoissonDrafting.h"
BOOST_AUTO_TEST_CASE(conditionalPoissonDrafting1, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonDraftingArgs args;
	std::vector<int>& indices = args.indices;
	std::vector<sampling::mpfr_class>& inclusionProbabilities = args.inclusionProbabilities;
	std::vector<sampling::mpfr_class>& weights = args.weights;

	boost::mt19937 randomSource;
	randomSource.seed(1);

	args.n = 2;
	int table[3];
	table[0] = table[1] = table[2] = 0;
	weights.push_back(3.0/6.0);
	weights.push_back(4.0/6.0);
	weights.push_back(5.0/6.0);
	double total = (3.0*4.0*1.0 + 3.0*2.0*5.0 + 3.0*4.0*5.0) / (6.0*6.0*6.0);
	double inclusion1 = (3.0*4.0*1.0 / (6.0*6.0*6.0)) + (3.0*2.0*5.0 / (6.0*6.0*6.0));
	double inclusion2 = (3.0*4.0*1.0 / (6.0*6.0*6.0)) + (3.0*4.0*5.0 / (6.0*6.0*6.0));
	double inclusion3 = (3.0*2.0*5.0 / (6.0*6.0*6.0)) + (3.0*4.0*5.0 / (6.0*6.0*6.0));
	double subset1 = (3.0*4.0*1.0 / (6.0*6.0*6.0)) / total;
	double subset2 = (3.0*2.0*5.0 / (6.0*6.0*6.0)) / total;
	double subset3 = (3.0*4.0*5.0 / (6.0*6.0*6.0)) / total;
	for(int i = 0; i < 5000; i++)
	{
		sampling::conditionalPoissonDrafting(args, randomSource);
		BOOST_TEST((int)indices.size() == 2);
		std::sort(indices.begin(), indices.end());
		if(indices[0] == 0 && indices[1] == 1)
		{
			BOOST_TEST(inclusionProbabilities[0].convert_to<double>() == inclusion1 / total);
			BOOST_TEST(inclusionProbabilities[1].convert_to<double>() == inclusion2 / total);
			table[0]++;
		}
		else if(indices[0] == 0 && indices[1] == 2)
		{
			BOOST_TEST(inclusionProbabilities[0].convert_to<double>() == inclusion1 / total);
			BOOST_TEST(inclusionProbabilities[2].convert_to<double>() == inclusion3 / total);
			table[1]++;
		}
		else if(indices[0] == 1 && indices[1] == 2)
		{
			BOOST_TEST(inclusionProbabilities[1].convert_to<double>() == inclusion2 / total);
			BOOST_TEST(inclusionProbabilities[2].convert_to<double>() == inclusion3 / total);
			table[2]++;
		}
	}
	BOOST_TEST(table[0]/5000.0 == subset1, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[1]/5000.0 == subset2, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[2]/5000.0 == subset3, boost::test_tools::tolerance(0.02));
}
BOOST_AUTO_TEST_CASE(conditionalPoissonDrafting2, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonDraftingArgs args;
	std::vector<int>& indices = args.indices;
	std::vector<sampling::mpfr_class>& weights = args.weights;

	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.push_back(2.0/5.0);
	weights.push_back(4.0/5.0);
	weights.push_back(4.0/5.0);
	args.n = 2;
	int table[3];
	table[0] = table[1] = table[2] = 0;
	long long size = 20000;
	for(int i = 0; i < size; i++)
	{
		conditionalPoissonDrafting(args, randomSource);
		std::sort(indices.begin(), indices.end());
		BOOST_TEST((int)indices.size() == 2);
		if(indices[0] == 0 && indices[1] == 1) table[0]++;
		if(indices[0] == 0 && indices[1] == 2) table[1]++;
		if(indices[0] == 1 && indices[1] == 2) table[2]++;
	}
	BOOST_TEST(table[0]/(double)size == 0.125, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[1]/(double)size == 0.125, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[2]/(double)size == 0.75, boost::test_tools::tolerance(0.02));
}
