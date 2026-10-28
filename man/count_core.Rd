\name{count_core}
\alias{count_core}
\title{Tuple counts for every set of t columns}
\usage{
count_core(x, t, vs = apply(x, 2, function(col) max(col, na.rm = TRUE) + 1L),
  rowmult = TRUE)
}
\arguments{
  \item{x}{an integer matrix (or a matrix or data frame that can be
    converted to one) with rows as runs and columns as factors, coded
    from 0. \code{NA} cells are allowed.}
  \item{t}{the number of columns per projection, an integer from 1 to
    \code{ncol(x)}.}
  \item{vs}{the number of levels of each column, a vector of length
    \code{ncol(x)}. The default is the largest symbol of each column
    plus one. Every symbol of column \code{j} must lie in
    \code{0:(vs[j]-1)}, otherwise an error is given.}
  \item{rowmult}{logical. If \code{TRUE} (the default), the matrix of
    row multiplicities is computed as well.}
}
\value{
  A list with elements
  \item{tots}{for every set of \code{t} columns, the number of possible
    value tuples, the product of their \code{vs}.}
  \item{ncovered}{for every set of \code{t} columns, the number of
    tuples that occur in at least one row.}
  \item{missing}{\code{tots - ncovered}.}
  \item{rowmult}{if \code{rowmult = TRUE}, an \code{nrow(x)} by
    \code{choose(ncol(x), t)} matrix whose entry for a row and a column
    set is the number of rows with the same tuple on that column set,
    or 0 if the row has an \code{NA} there. \code{NULL} otherwise.}
  The column sets are taken in lexicographic order, the order of
  \code{combn(ncol(x), t)}. The first three elements are numeric
  vectors of length \code{choose(ncol(x), t)}.
}
\description{
  The counting kernel behind coverage checks of covering arrays. For
  every set of \code{t} columns, in lexicographic order, the function
  counts how often each value tuple occurs among the rows. A row with
  an \code{NA} in a projection contributes nothing to that projection.
}
\details{
  The order of the column sets matches \code{combn} and the column
  sets reported by \code{CAs::coverage}. Tuples are linearized with
  the leftmost column varying fastest. The tuple space of one column
  set is limited to \eqn{2^{31}} values, and the \code{rowmult} matrix
  to \eqn{2 \times 10^9}{2e9} cells. The computation is single
  threaded and can be interrupted.
}
\seealso{
  \code{coverage} in package \pkg{CAs}.
}
\examples{
## the full factorial 2^3 covers every pair, each pair tuple twice
x <- as.matrix(expand.grid(0:1, 0:1, 0:1))
r <- count_core(x, 2)
r$ncovered
r$missing
r$rowmult[1:3, ]

## four of the runs form an orthogonal array of strength 2
oa <- x[c(1, 4, 6, 7), ]
count_core(oa, 2)$missing
## without its first run, one pair is missing in every two columns
count_core(oa[-1, ], 2, vs = c(2, 2, 2))$missing
}
