#include "setDefaultPrec.h"
#include "includeMPFRSampling.h"
namespace sampling
{
	SEXP setDefaultPrec(SEXP prec_sexp)
	{
	BEGIN_RCPP
		sampling::mpfr_class::default_precision(Rcpp::as<int>(prec_sexp));
	VOID_END_RCPP
		return R_NilValue;
	}
}
