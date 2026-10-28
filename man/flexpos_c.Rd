\name{flexpos_c}
\alias{flexpos_c}
\alias{markflex_c}
\title{Flexible positions and greedy marking of flexible values}
\usage{
flexpos_c(D, t, conv = c("flexpos", "coverage"))

markflex_c(D, t, fixrows = 0)
}
\arguments{
  \item{D}{a matrix or data frame with rows as runs and columns as
    factors, coded from 0 or from 1 with consecutive symbols in each
    column. \code{NA} marks a flexible value. The number of levels of
    each column is its largest symbol plus one (after a shift of codes
    from 1). Codes that are not consecutive from 0 or from 1 give an
    error.}
  \item{t}{the strength, an integer from 2 to \code{ncol(D)}.}
  \item{conv}{the reading of \code{NA} cells in \code{flexpos_c}.
    \code{"flexpos"} (the default) is the convention of
    \code{CAs::flexpos}: \code{NA} is an ordinary symbol that matches
    only \code{NA}. With \code{"coverage"} a row whose projection
    contains an \code{NA} contributes nothing to that projection and
    imposes no constraint on its cells. The two agree on arrays
    without \code{NA}, and \code{"flexpos"} never flags a cell that
    \code{"coverage"} would not flag.}
  \item{fixrows}{an integer from 0 to \code{nrow(D)}. The first
    \code{fixrows} rows are not moved.}
}
\value{
  \code{flexpos_c} returns a logical matrix of the dimensions of
  \code{D}, \code{TRUE} at flexible positions. Cells that are
  \code{NA} in \code{D} are \code{TRUE}.

  \code{markflex_c} returns a matrix in the coding of \code{D}, with
  the rows in their sorted order and the flexible values set to
  \code{NA}. It has the attributes \code{rowOrder} (the indices of the
  rows of \code{D} in the order returned) and \code{flexcount} (the
  number of \code{NA} cells in each returned row).
}
\description{
  C implementations of \code{flexpos} and \code{markflex} from the
  package \pkg{CAs}. A position of a covering array is flexible when
  every \code{t}-tuple through it is also covered by another row, so
  the position can be changed on its own without loss of coverage.
  \code{markflex_c} applies the greedy marking of Nayeri, Colbourn and
  Konjevod (2013), which turns a set of positions into \code{NA}
  (\dQuote{don't care}) values that are flexible together.
}
\details{
  A position counts as flexible in \code{flexpos_c} when, for every
  set of \code{t} columns through its column, the row's tuple on that
  set occurs in at least two rows.

  \code{markflex_c} sorts the rows by their number of flexible
  positions in ascending order (a stable sort, keeping the first
  \code{fixrows} rows in place, and no sort at all if
  \code{fixrows >= nrow(D) - 1}), as \code{CAs::markflex} does. It
  then goes through the column sets in lexicographic order and the
  rows from the top. The first row that covers a tuple keeps its
  \code{t} cells of that set, and a row with an \code{NA} in the
  projection covers nothing. Cells that are kept by no column set
  become \code{NA}.

  Both functions reproduce \code{CAs::flexpos} and
  \code{CAs::markflex} cell for cell on the test batteries of the
  package. The differences are in the input handling. Codes from 1 are
  shifted to 0 internally and back on output, codes that are not
  consecutive give an informative error, and \code{markflex_c} returns
  the attributes \code{rowOrder} and \code{flexcount}, where
  \code{CAs::markflex} returns \code{rowOrder} and \code{flexible}.
}
\references{
  Nayeri, P., Colbourn, C.J., Konjevod, G. (2013). Randomized
  post-optimization of covering arrays. \emph{European Journal of
  Combinatorics} \bold{34}, 91--103. \doi{10.1016/j.ejc.2012.07.017}
}
\seealso{
  \code{flexpos} and \code{markflex} in package \pkg{CAs},
  \code{\link{postopNCK_c}}.
}
\examples{
plan <- rbind(c(2,2,2,2,0,0,1), c(0,0,2,2,2,2,2), c(1,1,2,1,2,0,1),
              c(0,1,2,0,1,1,0), c(2,2,1,0,2,0,2), c(2,1,0,2,1,2,2),
              c(1,2,0,0,0,2,0), c(0,2,1,1,1,2,1), c(1,0,1,2,1,0,0),
              c(2,0,0,1,2,1,0), c(0,0,0,0,0,0,1), c(1,1,1,1,0,1,2),
              c(2,2,2,2,2,2,2), c(2,2,2,2,0,1,1))
## positions not involved in any pair that only one row covers
flexpos_c(plan, 2)
## a set of positions that can be flexible together
## the last two rows are entirely flexible and can be dropped
(m <- markflex_c(plan, 2))
attr(m, "flexcount")
}
