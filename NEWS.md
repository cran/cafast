# cafast 0.1.1

* New function `CA_to_PCA_c()`: brings a strength 2 array into the
  partitioned form used by `CAs::productPCA()`, choosing the `v` rows
  that give the widest left block (the most columns showing all `v`
  levels) by an exact search in compiled code. For the table's
  CA(25; 2, 220, 3) it gives `k1 = 137`, where `CAs::CA_to_PCA()` gives 64.

* `maxconstant_c()` counts a repeated constant row once. Constant rows
  already present are taken one per symbol; further copies stay in the
  array as ordinary rows and are never reported as constant rows (with
  `remove = TRUE` they are not removed). Before, duplicates were counted
  in two places: with `v` constant rows present they were returned
  without checking that they were `v` distinct constants (an array with
  the all-zero row twice and the all-one row once was reported as having
  3 constant rows at v = 3), and when the largest clique happened to
  have as many rows as the constant rows present, the duplicated rows
  were reported as the constant set (rows 000, 000, 012, 120, 211 at
  v = 3 gave rows 1 and 2 instead of a valid pair such as 1 and 5).
  The first was found 2026-10-20 by differential testing against the
  standalone `constant_rows` program in ca-tools, the second 2026-10-21
  by reading the code. `CAs::maxconstant()`
  has both shortcuts (checked with CAs 0.24-4).

# cafast 0.1.0

First release.

* `maxconstant_c()`, a C implementation of `CAs::maxconstant()`. Arrays
  with `NA` cells are accepted and treated as in the current
  `CAs::maxconstant()` (version 0.24-4): an `NA` cell matches no symbol
  when rows are compared, and the `NA` cells of the rows made constant
  are filled.
* `powerCT_any()`, `power_plan()` and `power_build()`, the power
  construction of Colbourn and Torres-Jimenez for any strength, number
  of columns and number of levels, homogeneous case.
* `count_core()`, tuple counts for every set of t columns.
* `flexpos_c()`, `markflex_c()` and `postopNCK_c()`, C versions of
  `CAs::flexpos()`, `CAs::markflex()` and of the marking steps of
  `CAs::postopNCK()`.
