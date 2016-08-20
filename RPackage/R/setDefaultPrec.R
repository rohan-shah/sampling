setDefaultPrec <- function(prec)
{
	.Call("setDefaultPrec", prec, PACKAGE="sampling")
}
