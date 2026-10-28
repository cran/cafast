\name{cafast-package}
\alias{cafast-package}
\alias{cafast}
\docType{package}
\title{Fast Companions to the CAs Package for Covering Arrays}
\description{
  Compiled implementations and generalizations of selected functions of
  the package \pkg{CAs} by Ulrike Groemping, for working with covering
  arrays.

  \itemize{
    \item \code{\link{maxconstant_c}}: makes the largest possible number
      of rows constant, as \code{CAs::maxconstant} does, with the clique
      search in C.
    \item \code{\link{powerCT_any}}, \code{\link{power_plan}},
      \code{\link{power_build}}: the power construction of Colbourn and
      Torres-Jimenez for any strength, number of columns and number of
      levels, where \code{CAs::powerCA} covers its catalogued settings.
    \item \code{\link{count_core}}: tuple counts for every set of
      \code{t} columns, the counting behind coverage checks.
    \item \code{\link{flexpos_c}}, \code{\link{markflex_c}}: flexible
      positions and greedy marking of flexible values, as
      \code{CAs::flexpos} and \code{CAs::markflex}.
    \item \code{\link{postopNCK_c}}: the post-optimization of Nayeri,
      Colbourn and Konjevod, as \code{CAs::postopNCK}.
  }

  Each function was tested against its counterpart in \pkg{CAs}. The
  help pages list the differences that remain, which concern
  tie-breaks, missing values and input checks. The test batteries,
  which download arrays from the NIST covering array tables, are
  installed in the \file{tests} folder of the installed package.
}
\author{
  Christopher Smolen \email{chrispbsmolen@gmail.com}
}
\seealso{
  The package \pkg{CAs}, \url{https://CRAN.R-project.org/package=CAs}.
}
\keyword{package}
