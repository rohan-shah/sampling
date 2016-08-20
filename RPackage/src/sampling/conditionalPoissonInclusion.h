#ifndef CONDITIONAL_POISSON_INCLUSION_RPACKAGE_HEADER_GUARD
#define CONDITIONAL_POISSON_INCLUSION_RPACKAGE_HEADER_GUARD
#include <Rcpp.h>
namespace sampling
{
	SEXP conditionalPoissonInclusion(SEXP sizes, SEXP n);
}
#endif
