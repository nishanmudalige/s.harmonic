\name{a.leg.poly}
\alias{a.leg.poly}
\title{A function to display the expression of an associated Legendre polynomials}
\description{Given a non-negative degree l and an order m, this function will
display the expression of the associated Legentre polynomial.
The output is an object of type 'expression'}
\usage{a.leg.poly(l, m)}
\arguments{
  \item{l}{The degree. Must be a non-negative integer.}
  \item{m}{The order. Must be an integer.}
}
\details{}
\value{
\item{}{The output is an object of tpye 'expression'}
}
\references{}
\author{Nishan Mudalige}
\note{}
\seealso{}
\examples{
## Display the expression of the associated Legendre polynomial of degree 3 and order 1
## a.leg.poly(3,1)

## Evaluate the associated Legendre polynomial of degree 5 and order -2
## a.leg.poly(5,-2)
}