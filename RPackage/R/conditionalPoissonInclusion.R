conditionalPoissonInclusion <- function(sizes, n)
{
	return(.Call("conditionalPoissonInclusion", sizes, n, PACKAGE="sampling"))
}
