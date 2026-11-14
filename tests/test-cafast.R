## Short checks of every exported function on small inline arrays.
## No network. The comparisons with 'CAs' run only if it is installed.
library(cafast)
check <- function(ok, what) if (!isTRUE(ok)) stop("test failed: ", what)
hasCAs <- requireNamespace("CAs", quietly = TRUE)

plan <- rbind(c(2,2,2,2,0,0,1), c(0,0,2,2,2,2,2), c(1,1,2,1,2,0,1),
              c(0,1,2,0,1,1,0), c(2,2,1,0,2,0,2), c(2,1,0,2,1,2,2),
              c(1,2,0,0,0,2,0), c(0,2,1,1,1,2,1), c(1,0,1,2,1,0,0),
              c(2,0,0,1,2,1,0), c(0,0,0,0,0,0,1), c(1,1,1,1,0,1,2),
              c(2,2,2,2,2,2,2), c(2,2,2,2,0,1,1))
nconst <- function(B) { # leading constant rows with ascending symbols
    cst <- !is.na(B[, 1]) & rowSums(B == B[, 1], na.rm = TRUE) == ncol(B); n <- 0L
    for (r in seq_len(nrow(B))) { if (!cst[r] || (r > 1 && B[r, 1] <= B[r - 1, 1])) break; n <- r }
    n }

## count_core: full factorial 2^3, every pair tuple twice
x <- as.matrix(expand.grid(0:1, 0:1, 0:1))
r <- count_core(x, 2)
check(identical(r$tots, c(4, 4, 4)) && identical(r$ncovered, c(4, 4, 4)) &&
      identical(r$missing, c(0, 0, 0)) && all(r$rowmult == 2), "count_core full factorial")
oa <- x[c(1, 4, 6, 7), ]   # OA(4, 3, 2, 2)
check(identical(count_core(oa, 2)$missing, c(0, 0, 0)) &&
      identical(count_core(oa[-1, ], 2, vs = c(2, 2, 2))$missing, c(1, 1, 1)),
      "count_core OA with one run removed")
xna <- x; xna[1, 1] <- NA
check(identical(count_core(xna, 2, vs = c(2, 2, 2))$rowmult[1, ], c(0, 0, 2)),
      "count_core NA row")

## maxconstant_c: rows 11 and 13 form the largest clique
M <- maxconstant_c(plan, verbose = 2)
check(nconst(M) == 2 && identical(attr(M, "constant_rows")$row_set_list, c(11L, 13L)),
      "maxconstant_c plan")
check(nrow(maxconstant_c(plan, remove = TRUE)) == 12, "maxconstant_c remove")
check(nconst(maxconstant_c(plan[-13, ], one_is_enough = TRUE)) == 1, "maxconstant_c one_is_enough")
## duplicated constant rows count once and stay in the array as ordinary rows
## (A: v = 2, the all-zero row twice; B: v = 3, the all-zero row twice and a
## largest clique of 2 that the duplicate pair would otherwise imitate)
A <- rbind(c(0,0,0), c(0,1,1), c(1,0,1), c(1,1,0), c(0,0,0))
B <- rbind(c(0,0,0), c(0,0,0), c(0,1,2), c(1,2,0), c(2,1,1))
MA <- maxconstant_c(A, verbose = 2); MB <- maxconstant_c(B, verbose = 2)
check(identical(attr(MA, "constant_rows")$row_set_list, 1L) && nconst(MA) == 1 && nrow(MA) == 5,
      "maxconstant_c duplicated constant row, v = 2")
check(identical(attr(MB, "constant_rows")$row_set_list, c(1L, 5L)) && nconst(MB) == 2 && nrow(MB) == 5,
      "maxconstant_c duplicated constant rows, v = 3")
check(nrow(maxconstant_c(B, remove = TRUE)) == 3, "maxconstant_c remove keeps the duplicate")
## OA(64, 4, 4, 3) with NA cells, one of them in a clique row
D <- as.matrix(expand.grid(0:3, 0:3, 0:3))
D <- unname(cbind(D, (D[, 1] + D[, 2] + D[, 3]) %% 4))
D[c(2, 22, 30), 4] <- NA
MD <- maxconstant_c(D, verbose = 2)
check(nconst(MD) == 4 && sum(is.na(MD)) == 2 && !anyNA(MD[1:4, ]), "maxconstant_c with NA")
check(identical(attr(MD, "constant_rows")$row_set_list, c(1L, 22L, 43L, 64L)), "maxconstant_c NA clique")
if (hasCAs) {
    check(identical(unname(maxconstant_c(plan)), unname(CAs::maxconstant(plan))), "maxconstant_c vs CAs, plan")
    ## NA cells in clique rows are handled by CAs from version 0.24-4 on
    if (utils::packageVersion("CAs") >= "0.24.4")
        check(nconst(CAs::maxconstant(D)) == nconst(MD), "maxconstant_c vs CAs, NA")
}

## flexpos_c and markflex_c
fp <- flexpos_c(plan, 2)
check(is.logical(fp) && all(dim(fp) == dim(plan)) && sum(fp) == 18, "flexpos_c count")
check(identical(fp, flexpos_c(plan, 2, conv = "coverage")), "flexpos_c conventions agree without NA")
check(identical(flexpos_c(plan + 1, 2), fp), "flexpos_c codes from 1")
mf <- markflex_c(plan, 2)
check(identical(as.numeric(attr(mf, "flexcount")), c(rep(0, 12), 7, 7)), "markflex_c flexcount")
check(identical(attr(mf, "rowOrder"), c(2L, 3L, 4L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 5L, 14L, 1L, 13L)),
      "markflex_c rowOrder")
if (hasCAs) {
    check(identical(unname(fp), unname(CAs::flexpos(plan, 2))), "flexpos_c vs CAs")
    h <- CAs::markflex(plan, 2)
    check(identical(as.numeric(mf), as.numeric(h)), "markflex_c vs CAs")
}

## postopNCK_c: two rows can go, coverage is kept
A <- postopNCK_c(plan, 2, seed = 1)
check(nrow(A) == 12 && attr(A, "seed") == 1, "postopNCK_c rows")
check(sum(count_core(A, 2, vs = rep(3, 7), rowmult = FALSE)$missing) == 0, "postopNCK_c coverage")

## power construction with a strength 2 binary ingredient from 'CAs'
if (hasCAs) {
    ing <- function(t, q, v) CAs::KSK(k = q)
    p <- power_plan(2, 100, 2, ingredient = ing, verbose = FALSE)
    check(all(p$N[1:3] == c(12, 12, 14)) && all(p$q[1:3] == c(11, 13, 5)), "power_plan")
    P <- powerCT_any(2, 100, 2, ingredient = ing, verbose = FALSE)
    check(identical(dim(P), c(12L, 100L)), "powerCT_any dim")
    check(sum(count_core(P, 2, vs = rep(2, 100), rowmult = FALSE)$missing) == 0, "powerCT_any coverage")
    B <- power_build(2, 100, 2, q = 11, e = 2, ingredient = CAs::KSK(k = 11))
    check(identical(dim(B), c(12L, 100L)) &&
          sum(count_core(B, 2, vs = rep(2, 100), rowmult = FALSE)$missing) == 0, "power_build")
}
## CA_to_PCA_c: the result is a PCA in the sense of CAs::is.PCA, still covers, and is at least as wide as CAs::CA_to_PCA
if (hasCAs) {
suppressMessages(library(CAs))
for (kv in list(c(10, 3), c(20, 4), c(12, 5))) {
    D <- unclass(CAs::bestCA(2, kv[1], kv[2])); D <- D - min(D)
    P <- CA_to_PCA_c(D); st <- attr(P, "PCAstatus")
    stopifnot(is.list(st), st$k1 + st$k2 == ncol(D), isTRUE(attr(P, "PCA_exact")))
    stopifnot(isTRUE(as.logical(CAs::is.PCA(P))), attr(CAs::is.PCA(P), "PCAstatus")$k1 == st$k1)
    stopifnot(all(CAs::coverage(P, 2)[c("total", "min")] == 1))
    old <- suppressMessages(attr(CAs::CA_to_PCA(D), "PCAstatus"))
    if (is.list(old)) stopifnot(st$k1 >= old$k1)
    ## brute force over all sets of v rows agrees on k1
    v <- kv[2]; if (choose(nrow(D), v) <= 20000) {
        b <- max(apply(utils::combn(nrow(D), v), 2, function(r) sum(apply(D[r, , drop = FALSE], 2, function(x) length(unique(x)) == v))))
        stopifnot(st$k1 == b)
    }
}
## the product of two results is a covering array of N + M - v rows
P1 <- CA_to_PCA_c(unclass(CAs::bestCA(2, 10, 3))); P2 <- CA_to_PCA_c(unclass(CAs::bestCA(2, 7, 3)))
Q <- CAs::productPCA(P1, P2)
check(nrow(Q) == nrow(P1) + nrow(P2) - 3 && all(CAs::coverage(Q, 2)[c("total", "min")] == 1), "CA_to_PCA_c product")
}

cat("cafast tests passed\n")
