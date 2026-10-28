\name{powerCT_any}
\alias{powerCT_any}
\alias{power_plan}
\alias{power_build}
\title{Power construction of covering arrays for any parameters}
\usage{
powerCT_any(t, k, v, emax = 4, qextra = 1, ingredient = NULL,
  verbose = TRUE, ...)

power_plan(t, k, v, emax = 4, qextra = 1, ingredient = NULL, maxq = 400,
  verbose = TRUE, ...)

power_build(t, k, v, q, e, ingredient = NULL, ...)
}
\arguments{
  \item{t}{the strength.}
  \item{k}{the number of columns.}
  \item{v}{the number of levels.}
  \item{emax}{the largest exponent \eqn{e} tried. Exponents
    \code{2:emax} are enumerated.}
  \item{qextra}{for each exponent, the number of admissible prime
    powers \eqn{q} tried beyond the smallest one.}
  \item{ingredient}{\code{NULL}, a matrix, or a function. \code{NULL}
    (the default) uses \code{CAs::bestCA(t, q, v, ...)}, which may
    download arrays from the NIST covering array tables. A function is
    called as \code{ingredient(t, q, v, ...)} and must return a
    CA(\eqn{t}, \eqn{q}, \eqn{v}) with at least \eqn{q} columns. A
    matrix is used as it is (only its first \eqn{q} columns), which is
    meant for \code{power_build}. \code{NA} cells of the ingredient are
    set to 0, and the ingredient is passed through
    \code{\link{maxconstant_c}}.}
  \item{maxq}{the largest prime power \eqn{q} considered.}
  \item{verbose}{logical. If \code{TRUE}, every admissible
    \eqn{(q, e)} and the resulting number of rows is printed.}
  \item{q}{a prime power, the number of levels of the orthogonal array
    and the number of columns of the ingredient.}
  \item{e}{the strength of the orthogonal array, at least 2.}
  \item{...}{passed to the ingredient function or to
    \code{CAs::bestCA}.}
}
\value{
  \code{power_plan} returns a data frame with one row per admissible
  \eqn{(q, e)}, sorted by the number of rows \code{N}, with columns
  \code{t}, \code{k}, \code{v}, \code{e}, \code{q}, \code{M} (rows of
  the hash family), \code{N1} (rows of the ingredient), \code{rho}
  (its constant rows), \code{chi} (constant rows added back),
  \code{N} and \code{kmax} (\eqn{q^e}, the largest number of columns
  of this \eqn{(q, e)}). The attribute \code{ingredients} is a list of
  the ingredient arrays by \eqn{q}. If no \eqn{(q, e)} is admissible,
  the result is \code{NULL}.

  \code{power_build} returns an integer matrix with \code{k} columns,
  levels \code{0:(v-1)} and no duplicated rows, with the attributes
  \code{t}, \code{q}, \code{e}, \code{M}, \code{N1} and \code{origin}.

  \code{powerCT_any} returns the array of \code{power_build} for the
  first row of the plan, with the plan as the further attribute
  \code{plan}. It gives a warning if the number of rows differs from
  the plan.
}
\description{
  The power construction of Colbourn and Torres-Jimenez (2010) for any
  strength \eqn{t}, number of columns \eqn{k} and number of levels
  \eqn{v}, in the homogeneous case with one ingredient. The function
  \code{powerCA} of the package \pkg{CAs} builds the 132 settings of
  its catalogue, and these functions extend the same construction to
  other settings.
}
\details{
  The columns of a Bush orthogonal array OA(\eqn{q^e}, \eqn{q+1},
  \eqn{q}, \eqn{e}) (from \code{CAs::SCA_Busht}), transposed, give a
  distributing hash family with \eqn{q^e} columns. Its first
  \eqn{M = (e-1)T(t, v) + 1} rows are used, where \eqn{T(t, v)} is the
  Turan number, so \eqn{q \ge M - 1} is needed, and \eqn{q^e \ge k}.
  An ingredient CA(\eqn{t}, \eqn{q}, \eqn{v}) with \eqn{N_1} rows,
  \eqn{\rho} of them constant after \code{\link{maxconstant_c}}, is
  placed in every row of the hash family, shifted by
  \eqn{(r-1)\rho} modulo \eqn{v} in row \eqn{r}. Then
  \eqn{\chi = \max(0, v - M(v - \rho))} constant rows are added back,
  which gives \eqn{N = M(N_1 - \rho) + \chi} rows. The assembly
  follows \code{CAs::DHHF2CA}.

  \code{power_plan} enumerates the admissible \eqn{(q, e)}, builds
  each ingredient once and reports \eqn{N}. \code{power_build} builds
  the array for one \eqn{(q, e)}. \code{powerCT_any} plans and builds
  the smallest. The constructions with several ingredients or with
  T-reductions in the catalogue of \code{CAs::powerCA} are not
  covered, and the result is often larger than the best known array.
}
\references{
  Colbourn, C.J., Torres-Jimenez, J. (2010). Heterogeneous hash
  families and covering arrays. In Bruen, A.A., Wehlau, D.L. (eds.),
  \emph{Contemporary Mathematics} \bold{523}, 3--15. American
  Mathematical Society, Providence, Rhode Island.
  \doi{10.1090/conm/523/10309}
}
\seealso{
  \code{powerCA}, \code{bestCA} and \code{DHHF2CA} in package
  \pkg{CAs}.
}
\examples{
## a strength 2 binary ingredient that needs no download
ing <- function(t, q, v) CAs::KSK(k = q)
power_plan(2, 100, 2, ingredient = ing)
A <- powerCT_any(2, 100, 2, ingredient = ing, verbose = FALSE)
dim(A)
sum(count_core(A, 2, vs = rep(2, 100), rowmult = FALSE)$missing)

B <- power_build(2, 100, 2, q = 11, e = 2, ingredient = CAs::KSK(k = 11))
dim(B)
}
