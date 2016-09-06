#include "conditionalPoissonInclusion.h"
#include "conditionalPoissonBase.h"
#include <vector>
namespace sampling
{
	SEXP conditionalPoissonInclusion(SEXP sizes_sexp, SEXP n_sexp)
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
		std::vector<double> inclusionProbabilities_double;
		std::transform(inclusionProbabilities.begin(), inclusionProbabilities.end(), std::back_inserter(inclusionProbabilities_double), [](mpfr_class x){ return x.convert_to<double>(); });
		return Rcpp::wrap(inclusionProbabilities_double);
	END_RCPP
	}
}
