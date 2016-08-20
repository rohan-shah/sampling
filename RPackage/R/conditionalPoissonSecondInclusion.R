conditionalPoissonSecondInclusion <- function(sizes, n)
{
	return(.Call("conditionalPoissonSecondInclusion", sizes, n, PACKAGE = "sampling"))
}
