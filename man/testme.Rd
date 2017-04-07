\name{testme}
\alias{testme}
\title{
Function for test the full capabilites of \code{DA.MRFA} package
}
\description{
testme is a function which includes all the functions included in the DA.MRFA package, and his purpose is to serve as a demonstration of the potential usage of the package functions.
}
\usage{
testme()
}
\details{
The function implements all the functions on the \code{DA.MRFA} package. For ilustrating purposes, the loading matrix obtained using MRFA (\code{\link{mrfa}} function) is rotated using the function \code{\link[PCovR]{promin}} (Lorenzo-Seva, 1999), from \code{PCovR} package (Vervloet, Kiers, Van den Noortgate & Ceulemans, 2015).
}
\references{
Lorenzo-Seva, U. (1999). Promin: A Method for Oblique Factor Rotation. Multivariate Behavioral Research, 34(3), 347–365 \url{http://doi.org/10.1207/S15327906MBR3403_3}

Vervloet, M., Kiers, H. A. L., Van den Noortgate, W., & Ceulemans, E. (2015). PCovR: An R Package for Principal Covariates Regression. JSS Journal of Statistical Software, 65(8). \url{http://doi.org/10.18637/jss.v065.i08}
}
\author{
David Navarro-Gonzalez

Urbano Lorenzo-Seva
}
\examples{
testme()
}
