\name{maxconstant_c}
\alias{maxconstant_c}
\title{Make the largest possible number of rows constant}
\usage{
maxconstant_c(D, verbose = 0, remove = FALSE, one_is_enough = FALSE,
  dupcheck = FALSE, ...)
}
\arguments{
  \item{D}{a numeric matrix, a uniform array with rows as runs and
    columns as factors. It should have no duplicated rows unless
    \code{dupcheck = TRUE}. \code{NA} cells are accepted (see Details),
    but a row consisting only of \code{NA} cells is an error.}
  \item{verbose}{one of 0, 1, 2 or 12, as in \code{CAs::maxconstant}.
    With 1 the rows used for the constant rows are reported as
    messages. With 2 the attribute \code{constant_rows} is added to the
    result. With 12 both happen. With 0 (the default) neither.}
  \item{remove}{logical. If \code{TRUE}, the constant rows are removed
    and the result has fewer rows than \code{D}. If \code{FALSE} (the
    default), they are moved to the front.}
  \item{one_is_enough}{logical. If \code{TRUE}, no clique search is
    made. Existing constant rows are moved to the front, and if there
    are none, the first row is made constant by swapping symbols within
    the columns.}
  \item{dupcheck}{logical. If \code{TRUE}, duplicated rows are removed
    before anything else is done.}
  \item{...}{currently not used.}
}
\value{
  A matrix equivalent to \code{D} (rows reordered, symbols permuted
  within columns) whose first rows are constant, with the levels of
  \code{D} in ascending order: the first constant row holds the
  smallest level, the second the next one, and so on. With
  \code{remove = TRUE} these rows are dropped instead. For
  \code{verbose} 2 or 12 the result carries the attribute
  \code{constant_rows}, a list with elements \code{design_name} (the
  deparsed argument) and \code{row_set_list} (the indices, in
  \code{D} after any \code{dupcheck}, of the rows that became the
  constant rows).
}
\description{
  A C implementation of \code{maxconstant} from the package \pkg{CAs}
  by Ulrike Groemping. A set of rows can be made constant by permuting
  symbols within columns exactly when every two of its rows differ in
  every column. The function finds a largest such set (a largest
  clique in the graph whose vertices are the rows and whose edges join
  rows that differ in every column), moves it to the front and relabels
  the symbols of each column so that these rows become constant.
}
\details{
  The graph is built and searched in C. The clique search is a
  branch and bound that stops as soon as it holds \eqn{v} rows, the
  number of levels, since no more rows can be constant. Existing
  constant rows are tried first as a forced start of the clique. If an
  unconstrained clique is strictly larger it is used, otherwise the
  clique containing the existing constant rows is kept, which is the
  preference of \code{CAs::maxconstant}. If the array already has
  \eqn{v} constant rows, or no larger clique exists, the existing
  constant rows are moved to the front and nothing else changes.

  Differences from \code{CAs::maxconstant} (version 0.24-4). When
  several largest cliques exist, \code{CAs::maxconstant} takes the last
  one returned by \code{igraph::largest_cliques} and this function
  takes the first one found in row order. The number of constant rows
  is the same. The arrays then differ by a row order and a symbol
  permutation, and on arrays with duplicated rows the number of rows
  outside the clique that happen to be constant as well can differ.
  The levels are taken from the whole array, whereas
  \code{CAs::maxconstant} takes them from the first column and stops
  unless every column has the same number of levels.

  Missing values are treated as in \code{CAs::maxconstant}. A row
  with an \code{NA} cell is not an existing constant row. Two rows are
  adjacent when they differ in every column in which both are
  non-missing, so an \code{NA} cell matches no symbol. After the
  clique rows are moved to the front, the \code{NA} cells among them
  are filled column by column, from the top, with the levels that do
  not occur among the clique rows in that column, in ascending order.
  The symbols are then permuted, and \code{NA} cells outside the clique
  rows stay \code{NA}. With \code{one_is_enough = TRUE} and no
  existing constant row, \code{NA} cells in the first row are left as
  they are, as in \code{CAs::maxconstant}.

  Two cases with \code{NA} cells give different results. First, the
  clique search here stops at \eqn{v} rows, while
  \code{igraph::largest_cliques} can return a clique of more than
  \eqn{v} rows when \code{NA} cells make extra rows compatible. No
  filling makes more than \eqn{v} rows constant, and
  \code{CAs::maxconstant} then returns fewer constant rows (often
  none), while this function returns \eqn{v}. Second, when only one
  row can be made constant, none is constant yet and the first row has
  an \code{NA} cell, \code{CAs::maxconstant} stops with an error and
  this function fills the cell. The fill uses the array's own levels,
  where \code{CAs::maxconstant} uses \code{0:(v-1)}, which is the same
  for arrays coded from 0.
}
\references{
  Torres-Jimenez, J., Acevedo-Juarez, B., Avila-George, H. (2021).
  Covering array EXtender. \emph{Applied Mathematics and Computation}
  \bold{402}, 126122. \doi{10.1016/j.amc.2021.126122}
}
\seealso{
  \code{maxconstant} in package \pkg{CAs}, \code{\link{powerCT_any}},
  which uses this function on its ingredient.
}
\examples{
## a 3-level array: one row is constant and one more can be made constant
plan <- rbind(c(2,2,2,2,0,0,1), c(0,0,2,2,2,2,2), c(1,1,2,1,2,0,1),
              c(0,1,2,0,1,1,0), c(2,2,1,0,2,0,2), c(2,1,0,2,1,2,2),
              c(1,2,0,0,0,2,0), c(0,2,1,1,1,2,1), c(1,0,1,2,1,0,0),
              c(2,0,0,1,2,1,0), c(0,0,0,0,0,0,1), c(1,1,1,1,0,1,2),
              c(2,2,2,2,2,2,2), c(2,2,2,2,0,1,1))
maxconstant_c(plan, verbose = 1)

## an orthogonal array OA(64, 4, 4, 3) with three NA cells, one of
## them in row 22, which becomes the second constant row
D <- as.matrix(expand.grid(0:3, 0:3, 0:3))
D <- unname(cbind(D, (D[, 1] + D[, 2] + D[, 3]) \%\% 4))
D[c(2, 22, 30), 4] <- NA
head(maxconstant_c(D, verbose = 1))
}
