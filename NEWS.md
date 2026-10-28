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
