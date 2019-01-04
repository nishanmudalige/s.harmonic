\name{spherical.harmonic}
\alias{spherical.harmonic}
\title{A function to calculate the real spherical harmonics}
\description{A function which calculates the real (tesseral) spherical harmonics for a given degree l and order m. The spherical harmonic is evaluated at the spherical coordinate (theta,phi) on the unit sphere in S^2.}
\usage{spherical.harmonic(l, m, theta, phi)}
\arguments{
  \item{l}{The degree of the associated Legendre polynomial. Must me a non negative integer.}
  \item{m}{The order of the associated Legendre polynomial. Must be an integer.}
  \item{theta}{The polar angle.}
  \item{phi}{The azimuthal angle.}
}
\details{}
\value{}
\references{}
\author{Nishan Mudalige}
\note{}
\seealso{}
\examples{
## spherical.harmonic(3, 1, pi/2, pi/4)

## spherical.harmonic(5, -2, pi/3, pi/6)
}
