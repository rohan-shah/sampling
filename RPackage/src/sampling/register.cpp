#include <Rcpp.h>
#include <internal.h>
#include "conditionalPoissonInclusion.h"
#include "conditionalPoissonSecondInclusion.h"
#include "setDefaultPrec.h"
#ifdef _MSC_VER
	#undef RcppExport
	#define RcppExport extern "C" __declspec(dllexport)
#endif
extern "C" const char* package_name = "sampling";
R_CallMethodDef callMethods[] = 
{
	{"conditionalPoissonInclusion", (DL_FUNC)&sampling::conditionalPoissonInclusion, 2},
	{"conditionalPoissonSecondInclusion", (DL_FUNC)&sampling::conditionalPoissonSecondInclusion, 2},
	{"setDefaultPrec", (DL_FUNC)&sampling::setDefaultPrec, 1},
	{NULL, NULL, 0}
};
RcppExport void R_init_sampling(DllInfo *info)
{
	std::vector<R_CallMethodDef> callMethodsVector;
	R_CallMethodDef* packageCallMethods = callMethods;
	while(packageCallMethods->name != NULL) packageCallMethods++;
	callMethodsVector.insert(callMethodsVector.begin(), callMethods, packageCallMethods);

	R_CallMethodDef* RcppStartCallMethods = Rcpp_get_call();
	R_CallMethodDef* RcppEndCallMethods = RcppStartCallMethods;
	while(RcppEndCallMethods->name != NULL) RcppEndCallMethods++;
	callMethodsVector.insert(callMethodsVector.end(), RcppStartCallMethods, RcppEndCallMethods);
	R_CallMethodDef blank = {NULL, NULL, 0};
	callMethodsVector.push_back(blank);

	R_registerRoutines(info, NULL, &(callMethodsVector[0]), NULL, NULL);
	init_Rcpp_cache();
}
