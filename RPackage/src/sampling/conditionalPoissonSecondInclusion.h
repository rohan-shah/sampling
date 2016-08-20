#ifndef CONDITIONAL_POISSON_SECOND_INCLUSION_RPACKAGE_HEADER_GUARD
#define CONDITIONAL_POISSON_SECOND_INCLUSION_RPACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace sampling
{
	SEXP conditionalPoissonSecondInclusion(SEXP sizes_sexp, SEXP n_sexp);
}
#endif
