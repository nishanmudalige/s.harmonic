\name{a.leg.poly}
\alias{a.leg.poly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function to express associated Legendre polynomials
%%  ~~function to do ... ~~
}
\description{Given a non-negative degree l and an order m, this function will
evaluate the expression of the associated Legentre polynomial.
The output is an object of type 'expression'
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
a.leg.poly(l, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{l}{The degree. Must be a non-negative integer.
%%     ~~Describe \code{x} here~~
}
  \item{m}{The order. Must be an integer.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{}{The output is an object of tpye 'expression'}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Nishan Mudalige
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
a.leg.poly(3,1)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.