\name{CA_to_PCA_c}
\alias{CA_to_PCA_c}
\title{Bring a strength 2 array into partitioned form with the widest left block}
\usage{
CA_to_PCA_c(D, seconds = 60, ...)
}
\arguments{
  \item{D}{a numeric matrix, a uniform covering array of strength 2 with
    rows as runs and columns as factors.}
  \item{seconds}{time limit for the search for the best rows. When it is
    reached, the best rows found so far are used and the attribute
    \code{PCA_exact} is \code{FALSE}.}
  \item{...}{currently not used.}
}
\value{
  A matrix equivalent to \code{D} (rows and columns reordered, symbols
  permuted within columns) in the partitioned form that
  \code{CAs::productPCA} uses: the first \code{v} rows read as the
  levels in ascending order in each of the first \code{k1} columns, and
  in the other \code{k2} columns the level in row \code{i} is not larger
  than the \code{i}th level. Attributes of \code{D} are kept. Added
  attributes: \code{PCAstatus}, a list with \code{type} (\code{"PCA"},
  or \code{"SCA"} when the top right block is constant at the smallest
  level or empty), \code{k1} and \code{k2}, as set by
  \code{CAs::CA_to_PCA}, or \code{FALSE} when no column can be in the
  left block; \code{PCA_rows}, the rows of \code{D} that became the top
  block; \code{PCA_columns}, the columns of \code{D} in their new order;
  and \code{PCA_exact}.
}
\description{
  Chooses the \code{v} rows of \code{D} on which the largest number of
  columns show all \code{v} levels, moves them to the top and those
  columns to the left, and permutes the symbols within every column.
  The number of columns of a product of two partitioned arrays grows
  with the width \code{k1} of this left block.
}
\details{
  \code{CAs::CA_to_PCA} works on the first \code{v} rows after
  \code{maxconstant}, or with \code{tryhard = TRUE} looks for a split
  with \code{k2} equal to 1 or 2. This function searches all sets of
  \code{v} rows instead, by depth-first search with a bound on the
  columns still usable, in compiled code. For three levels and a few
  hundred columns it takes seconds. The result is exact unless the time
  limit is reached.

  In the columns of the right block the symbols are numbered in the
  order in which they first appear from the top, which gives the
  condition on row \code{i} above.
}
\references{
  Colbourn, C. J., Martirosyan, S. S., Mullen, G. L., Shasha, D.,
  Sherwood, G. B. and Yucas, J. L. (2006). Products of mixed covering
  arrays of strength two. \emph{Journal of Combinatorial Designs}
  \bold{14}, 124--138. \doi{10.1002/jcd.20065}
}
\seealso{\code{CAs::CA_to_PCA}, \code{CAs::productPCA}, \code{CAs::is.PCA}}
\examples{
## the orthogonal array OA(9; 2, 4, 3) from two coordinates and two sums
D <- as.matrix(expand.grid(0:2, 0:2))
D <- unname(cbind(D, (D[, 1] + D[, 2]) \%\% 3, (D[, 1] + 2 * D[, 2]) \%\% 3))
P <- CA_to_PCA_c(D)
attr(P, "PCAstatus")
P[1:3, ]
}
