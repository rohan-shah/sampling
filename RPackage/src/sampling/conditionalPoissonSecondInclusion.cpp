#include "conditionalPoissonSecondInclusion.h"
#include "conditionalPoissonBase.h"
namespace sampling
{
	SEXP conditionalPoissonSecondInclusion(SEXP sizes_sexp, SEXP n_sexp)
	{
	BEGIN_RCPP
		Rcpp::NumericVector sizes = Rcpp::as<Rcpp::NumericVector>(sizes_sexp);
		int n = Rcpp::as<int>(n_sexp);
		conditionalPoissonArgs args;
		args.weights.insert(args.weights.begin(), sizes.begin(), sizes.end());
		args.n = n;

		std::vector<mpfr_class> inclusionProbabilities;
		args.zeroWeights.resize(sizes.size());
		args.deterministicInclusion.resize(sizes.size());
		args.indices.clear();
		for(int i = 0; i != sizes.size(); i++)
		{
			args.zeroWeights[i] = args.deterministicInclusion[i] = false;
			if(sizes[i] < 0 || sizes[i] > 1)
			{
				throw std::runtime_error("Sizes must be values in [0, 1]");
			}
			else if(sizes[i] == 0)
			{
				args.zeroWeights[i] = true;
			}
			else if(sizes[i] == 1)
			{
				args.deterministicInclusion[i] = true;
				args.indices.push_back(i);
			}
		}
		computeExponentialParameters(args);
		conditionalPoissonInclusionProbabilities(args, inclusionProbabilities);
		boost::numeric::ublas::matrix<mpfr_class> secondOrder(sizes.size(), sizes.size());
		conditionalPoissonSecondOrderInclusionProbabilities(args, inclusionProbabilities, secondOrder);
		Rcpp::NumericMatrix result(sizes.size(), sizes.size());
		for(int i = 0; i < sizes.size(); i++)
		{
			for(int j = 0; j < sizes.size(); j++)
			{
				result(i, j) = secondOrder(i, j).convert_to<double>();
			}
		}
		return result;
	END_RCPP
	}
}
