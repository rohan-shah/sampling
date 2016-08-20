.onLoad <- function(libname, pkgname)
{
	library.dynam(package="sampling", chname="sampling", lib.loc = .libPaths())
}
