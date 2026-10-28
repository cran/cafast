\name{postopNCK_c}
\alias{postopNCK_c}
\title{Randomized post-optimization of a covering array}
\usage{
postopNCK_c(D, t, fixrows = 0, outerRetry = 50, outerMaxnochange = 10,
  innerRetry = 10, innerMaxnochange = 25, seed = NULL,
  timeBudget = Inf, restarts = 0, verbose = 0)
}
\arguments{
  \item{D}{a covering array of strength \code{t}, a matrix or data
    frame with rows as runs, coded from 0 or from 1 with consecutive
    symbols in each column. Every column must contain all of its
    levels. \code{NA} marks a flexible value.}
  \item{t}{the strength, an integer from 2 to \code{ncol(D)}.}
  \item{fixrows}{an integer from 0 to \code{nrow(D)}. The first
    \code{fixrows} rows are never moved or removed.}
  \item{outerRetry}{the largest number of outer rounds.}
  \item{outerMaxnochange}{the number of consecutive outer rounds that
    remove no row after which the search stops (or restarts, see
    \code{restarts}).}
  \item{innerRetry}{the number of inner retries in one outer round.
    Each retry works on a different last row.}
  \item{innerMaxnochange}{the number of inner iterations without an
    increase in the number of flexible values of the last row after
    which the next inner retry starts.}
  \item{seed}{\code{NULL} or an integer seed for the random process. If
    \code{NULL}, a seed is drawn with \code{sample(32000, 1)}. The seed
    is passed to \code{set.seed} and stored in the attribute
    \code{seed} of the result.}
  \item{timeBudget}{a time limit in seconds (default \code{Inf}). When
    it expires, the best array found so far is returned with the
    attribute \code{timeout = TRUE}.}
  \item{restarts}{the number of restarts allowed after the stopping
    rule of \code{outerMaxnochange} is reached. A restart fills every
    \code{NA} of the best array at random, marks it again and
    continues. The default 0 keeps the stopping behaviour of
    \code{CAs::postopNCK}.}
  \item{verbose}{0 (the default) prints nothing. With 1 or more, the
    number of rows is printed in every outer round, and with 2 or more
    every inner retry is printed as well.}
}
\value{
  The array with the fewest rows found, of class \code{"ca"} added to
  the matrix classes, in the coding of \code{D}, with its remaining
  flexible values as \code{NA}. It has the attribute \code{seed}, the
  attribute \code{flexible} (\code{list(value = NA)}) if it contains
  \code{NA} cells, and \code{timeout = TRUE} or \code{interrupt = TRUE}
  if the search was stopped by \code{timeBudget} or by a user
  interrupt. If \code{nrow(D)} already equals the product of the
  \code{t} largest numbers of levels, \code{D} is returned unchanged.
}
\description{
  The post-optimization of Nayeri, Colbourn and Konjevod (2013) as
  implemented in \code{postopNCK} of the package \pkg{CAs}, with the
  marking steps in C (\code{\link{flexpos_c}},
  \code{\link{markflex_c}}). Rows whose cells can all be made flexible
  are removed. The other flexible values are filled at random and the
  array is marked again, repeatedly, so that more rows become
  removable.
}
\details{
  A row with fewer than \code{t} non-flexible cells covers nothing and
  is removed. The inner loop moves a row with the most flexible values
  to the bottom, fills the flexible values (in a column where the
  bottom row has a symbol, with that symbol, elsewhere at random),
  marks the array again and repeats. Inner retries and outer rounds
  follow \code{CAs::postopNCK}.

  The function was calibrated against \code{CAs::postopNCK} by the
  distribution of the achieved number of rows over many seeds, not by
  trajectory. The same seed gives different runs in the two functions.
  The deliberate differences, agreed with the author of \pkg{CAs}, are
  these. Random fills use the array's own coding. Exhausted rows are
  kept in an explicit vector at the top of the array, and exactly
  these rows are frozen between inner retries (\code{CAs::postopNCK}
  freezes up to the last exhausted index, which in practice freezes
  the rows with flexible values as well). An inner pass that removes
  nothing hands its perturbed array to the next pass instead of
  returning its input. The function is silent by default. The time
  budget and the restarts are additions taken from the paper. An
  interrupt returns the best array found so far.
}
\references{
  Nayeri, P., Colbourn, C.J., Konjevod, G. (2013). Randomized
  post-optimization of covering arrays. \emph{European Journal of
  Combinatorics} \bold{34}, 91--103. \doi{10.1016/j.ejc.2012.07.017}
}
\seealso{
  \code{postopNCK} in package \pkg{CAs}, \code{\link{markflex_c}}.
}
\examples{
plan <- rbind(c(2,2,2,2,0,0,1), c(0,0,2,2,2,2,2), c(1,1,2,1,2,0,1),
              c(0,1,2,0,1,1,0), c(2,2,1,0,2,0,2), c(2,1,0,2,1,2,2),
              c(1,2,0,0,0,2,0), c(0,2,1,1,1,2,1), c(1,0,1,2,1,0,0),
              c(2,0,0,1,2,1,0), c(0,0,0,0,0,0,1), c(1,1,1,1,0,1,2),
              c(2,2,2,2,2,2,2), c(2,2,2,2,0,1,1))
A <- postopNCK_c(plan, 2, seed = 1)
dim(A)
## still covering: no pair is missing in any two columns
sum(count_core(A, 2, vs = rep(3, 7), rowmult = FALSE)$missing)
}
