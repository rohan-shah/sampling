#include <boost/test/unit_test.hpp>
#include "conditionalPoissonSequential.h"
BOOST_AUTO_TEST_CASE(conditionalPoissonSequential1, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonSequentialArgs args(false);
	std::vector<int>& indices = args.indices;
	std::vector<sampling::mpfr_class>& weights = args.weights;
	std::vector<sampling::mpfr_class>& inclusionProbabilities = args.inclusionProbabilities;
	
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
	long long size = 50000;
	for(int i = 0; i < size; i++)
	{
		sampling::conditionalPoissonSequential(args, randomSource);
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
	BOOST_TEST(table[0]/(double)size == subset1, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[1]/(double)size == subset2, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[2]/(double)size == subset3, boost::test_tools::tolerance(0.02));
}
BOOST_AUTO_TEST_CASE(conditionalPoissonSequential2, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonSequentialArgs args(false);
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
BOOST_AUTO_TEST_CASE(conditionalPoissonSequential3, * boost::unit_test::tolerance(0.00001))
{
	sampling::conditionalPoissonSequentialArgs args(false);
	std::vector<int>& indices = args.indices;
	std::vector<sampling::mpfr_class>& weights = args.weights;
	std::vector<sampling::mpfr_class>& inclusionProbabilities = args.inclusionProbabilities;
	
	boost::mt19937 randomSource;
	randomSource.seed(1);

	weights.push_back(0);
	weights.push_back(3.0/6.0);
	weights.push_back(0);
	weights.push_back(4.0/6.0);
	weights.push_back(0);
	weights.push_back(5.0/6.0);
	args.n = 2;
	int table[3];
	table[0] = table[1] = table[2] = 0;
	double total = (3.0*4.0*1.0 + 3.0*2.0*5.0 + 3.0*4.0*5.0) / (6.0*6.0*6.0);
	double inclusion1 = (3.0*4.0*1.0 / (6.0*6.0*6.0)) + (3.0*2.0*5.0 / (6.0*6.0*6.0));
	double inclusion2 = (3.0*4.0*1.0 / (6.0*6.0*6.0)) + (3.0*4.0*5.0 / (6.0*6.0*6.0));
	double inclusion3 = (3.0*2.0*5.0 / (6.0*6.0*6.0)) + (3.0*4.0*5.0 / (6.0*6.0*6.0));
	double subset1 = (3.0*4.0*1.0 / (6.0*6.0*6.0)) / total;
	double subset2 = (3.0*2.0*5.0 / (6.0*6.0*6.0)) / total;
	double subset3 = (3.0*4.0*5.0 / (6.0*6.0*6.0)) / total;
	long long size = 50000;
	for(int i = 0; i < size; i++)
	{
		conditionalPoissonSequential(args, randomSource);
		BOOST_TEST((int)indices.size() == 2);
		std::sort(indices.begin(), indices.end());
		if(indices[0] == 1 && indices[1] == 3)
		{
			BOOST_TEST(inclusionProbabilities[1].convert_to<double>() == inclusion1 / total);
			BOOST_TEST(inclusionProbabilities[3].convert_to<double>() == inclusion2 / total);
			table[0]++;
		}
		else if(indices[0] == 1 && indices[1] == 5)
		{
			BOOST_TEST(inclusionProbabilities[1].convert_to<double>() == inclusion1 / total);
			BOOST_TEST(inclusionProbabilities[5].convert_to<double>() == inclusion3 / total);
			table[1]++;
		}
		else if(indices[0] == 3 && indices[1] == 5)
		{
			BOOST_TEST(inclusionProbabilities[3].convert_to<double>() == inclusion2 / total);
			BOOST_TEST(inclusionProbabilities[5].convert_to<double>() == inclusion3 / total);
			table[2]++;
		}
		BOOST_TEST((int)indices.size() == 2);
	}
	BOOST_TEST(table[0]/(double)size == subset1, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[1]/(double)size == subset2, boost::test_tools::tolerance(0.02));
	BOOST_TEST(table[2]/(double)size == subset3, boost::test_tools::tolerance(0.02));
}

