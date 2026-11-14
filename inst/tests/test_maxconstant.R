suppressMessages({library(CAs); library(caverify)})
suppressMessages(library(cafast))
nconst <- function(B) { # leading constant rows with strictly ascending symbol = the clique placed in front (a chance duplicate behind it does not count)
    cst <- !is.na(B[, 1]) & rowSums(B == B[, 1], na.rm = TRUE) == ncol(B); n <- 0L
    for (r in seq_len(nrow(B))) { if (!cst[r] || (r > 1 && B[r, 1] <= B[r - 1, 1])) break; n <- r }
    n }
ntotal <- function(B) sum(!is.na(B[, 1]) & rowSums(B == B[, 1], na.rm = TRUE) == ncol(B))
mk <- function(t, k, v) { A <- suppressWarnings(nistCA(t, k, v)); A <- as.matrix(A); A[is.na(A)] <- 0L; storage.mode(A) <- "integer"; A }
specs <- list(c(2,20,4), c(3,30,3), c(2,40,5), c(3,20,4), c(4,15,3), c(3,60,4), c(2,60,6), c(4,30,4), c(5,20,3), c(4,20,5))
cases <- lapply(specs, function(s) mk(s[1], s[2], s[3]))
names(cases) <- sapply(specs, function(s) sprintf("nistCA(%d,%d,%d)", s[1], s[2], s[3]))
set.seed(7)
for (i in 1:40) { N <- sample(10:80, 1); k <- sample(3:12, 1); v <- sample(2:6, 1)
    M <- matrix(sample(0:(v - 1), N * k, TRUE), N, k); M[1:v, ] <- matrix(rep(0:(v - 1), k), v, k); cases[[sprintf("random%d", i)]] <- M[sample(N), ] }
cases[["bush64"]] <- { D <- lhs::createBush(4, 5); storage.mode(D) <- "integer"; D }
cases[["cyc19"]] <- { D <- cyc(19, 2); D <- as.matrix(D); storage.mode(D) <- "integer"; D }
fails <- 0; tot_h <- 0; tot_c <- 0
for (nm in names(cases)) {
    A <- cases[[nm]]; t_use <- if (grepl("nistCA", nm)) as.integer(sub("nistCA\\((\\d).*", "\\1", nm)) else 2L
    th <- system.time(H <- maxconstant(A))[["elapsed"]]
    tc <- system.time(C <- maxconstant_c(A))[["elapsed"]]
    tot_h <- tot_h + th; tot_c <- tot_c + tc
    same_n <- nconst(H) == nconst(C)
    front <- nconst(C) >= 1
    asc <- nconst(C) < 2 || all(diff(C[seq_len(nconst(C)), 1]) > 0)  # holds by definition of nconst; kept as a printed reminder
    equiv <- nrow(C) == nrow(A) && all(sapply(seq_len(ncol(A)), function(j) identical(sort(as.vector(table(A[, j]))), sort(as.vector(table(C[, j]))))))
    cov <- if (grepl("nistCA|bush|cyc", nm)) { r <- ca_verify(C, t_use); isTRUE(r$covered) } else NA
    ok <- same_n && front && asc && equiv && !isFALSE(cov)
    if (!ok) fails <- fails + 1
    cat(sprintf("%-16s %5d x %2d  hers %d front (%d total) %6.2fs   ours %d front (%d total) %6.2fs   asc %s equiv %s covers %s  %s\n",
        nm, nrow(A), ncol(A), nconst(H), ntotal(H), th, nconst(C), ntotal(C), tc, asc, equiv, cov, if (ok) "OK" else "MISMATCH"))
}
for (nm in c("bush64", "cyc19")) { A <- cases[[nm]]; C1 <- maxconstant_c(A, remove = TRUE); C2 <- maxconstant_c(A, one_is_enough = TRUE)
    cat(sprintf("%s: remove=TRUE gives %d rows (was %d); one_is_enough gives %d constant rows\n", nm, nrow(C1), nrow(A), nconst(C2))) }

## NA-bearing cases (added 2026-10-17, against her fixed maxconstant).
## blank(): set about frac of the cells to NA, keeping every level in every
## column (her symmetry check counts non-missing levels) and no row all NA.
blank <- function(A, frac, rows = seq_len(nrow(A)), seed) {
    set.seed(seed); v <- length(unique(as.vector(A)))
    cells <- which(row(A) %in% rows); nb <- max(1L, round(frac * length(cells)))
    for (cell in sample(cells)) {
        if (nb == 0) break
        B <- A; B[cell] <- NA; j <- (cell - 1L) %/% nrow(A) + 1L; r <- (cell - 1L) %% nrow(A) + 1L
        if (length(unique(na.omit(B[, j]))) < v || all(is.na(B[r, ]))) next
        A <- B; nb <- nb - 1L
    }
    A }
## check that C is A with rows c(cl, rest), each column relabelled by a
## bijection on the non-missing cells, NA kept outside the clique rows and
## filled inside them (cl = attr(C, "constant_rows")$row_set_list)
equivNA <- function(A, C, cl) {
    if (nrow(C) != nrow(A)) return(FALSE)
    A2 <- A[c(cl, setdiff(seq_len(nrow(A)), cl)), , drop = FALSE]; top <- seq_along(cl)
    if (anyNA(C[top, ])) return(FALSE)
    if (!identical(is.na(A2[-top, , drop = FALSE]), is.na(C[-top, , drop = FALSE]))) return(FALSE)
    all(sapply(seq_len(ncol(A)), function(j) { ok <- !is.na(A2[, j])
        m <- unique(cbind(A2[ok, j], C[ok, j])); !anyDuplicated(m[, 1]) && !anyDuplicated(m[, 2]) })) }
nacases <- list()
for (i in 1:40) { nm <- sprintf("random%d", i); nacases[[paste0(nm, "+NA")]] <- blank(cases[[nm]], 0.04, seed = 100 + i) }
## random arrays without planted constant rows: each level placed in
## random rows per column, so no clique is planted and many maxima are below v
set.seed(17)
for (i in 1:20) { N <- sample(10:60, 1); k <- sample(3:10, 1); v <- sample(3:6, 1)
    M <- matrix(sample(0:(v - 1), N * k, TRUE), N, k)
    for (j in seq_len(k)) M[sample(N, v), j] <- 0:(v - 1)
    nacases[[sprintf("free%d+NA", i)]] <- blank(M, 0.04, seed = 200 + i) }
## NAs inside the clique rows: blank cells of the planted constant rows of
## the random cases (they stop being constant, so her clique search runs
## and the NA cells of the clique must be filled), and of Bush's clique
for (i in 1:20) { nm <- sprintf("random%d", i); A <- cases[[nm]]
    cst <- which(rowSums(A == A[, 1]) == ncol(A))
    nacases[[paste0(nm, "+NAclique")]] <- blank(A, 0.25, rows = cst, seed = 300 + i) }
B <- cases[["bush64"]]; clB <- attr(maxconstant_c(B, verbose = 2), "constant_rows")$row_set_list
nacases[["bush64+NAclique"]] <- blank(B, 0.2, rows = clB, seed = 401)
nacases[["bush64+NA"]] <- blank(B, 0.03, seed = 402)
nafails <- 0; nabig <- 0; nadup <- 0
for (nm in names(nacases)) {
    A <- nacases[[nm]]
    th <- system.time(H <- maxconstant(A, verbose = 2))[["elapsed"]]
    tc <- system.time(C <- maxconstant_c(A, verbose = 2))[["elapsed"]]
    tot_h <- tot_h + th; tot_c <- tot_c + tc
    same_n <- nconst(H) == nconst(C)
    eqC <- equivNA(A, C, attr(C, "constant_rows")$row_set_list)
    eqH <- equivNA(A, H, attr(H, "constant_rows")$row_set_list)
    ## known difference: NA cells can make her igraph clique larger than v
    ## (v + 1 rows pairwise compatible); no filling makes more than v rows
    ## constant, hers then returns fewer constant rows, ours stops at v
    v <- length(unique(na.omit(as.vector(A))))
    big <- length(attr(H, "constant_rows")$row_set_list) > v
    ## known difference (cafast 0.1.1): hers can report copies of one constant
    ## row as several constant rows; ours counts a repeated constant row once,
    ## so it must give at least as many distinct constant rows as hers
    rsH <- attr(H, "constant_rows")$row_set_list
    dupH <- !big && length(rsH) > 1 && all(!is.na(A[rsH, 1])) &&
        all(rowSums(A[rsH, , drop = FALSE] == A[rsH, 1], na.rm = TRUE) == ncol(A)) && anyDuplicated(A[rsH, 1]) > 0
    ok <- if (big) nconst(C) == v && eqC else if (dupH) nconst(C) >= length(unique(A[rsH, 1])) && eqC else same_n && nconst(C) >= 1 && eqC
    if (big) nabig <- nabig + 1
    if (dupH) nadup <- nadup + 1
    if (!ok) nafails <- nafails + 1
    cat(sprintf("%-18s %3d x %2d  NA %3d (clique rows %d)  hers %d front  ours %d front  equiv ours %s hers %s  %s\n",
        nm, nrow(A), ncol(A), sum(is.na(A)), sum(is.na(A[attr(C, "constant_rows")$row_set_list, ])),
        nconst(H), nconst(C), eqC, eqH, if (!ok) "MISMATCH" else if (big) "OK (her clique > v)" else if (dupH) "OK (her rows repeat a constant)" else "OK"))
}
fails <- fails + nafails
## documented differences, not counted: hers fails, ours gives a usable array
D1 <- as.matrix(cyc(19, 2)); D1[1, 3] <- NA   # only one constant row possible, NA in it
r1 <- tryCatch({ maxconstant(D1); "no error" }, error = function(e) conditionMessage(e))
cat(sprintf("single-row clique with NA: hers %s; ours %d constant row\n", r1, nconst(maxconstant_c(D1))))
D2 <- rbind(c(0, NA), c(1, NA), c(NA, 0), c(NA, 1), c(0, 1), c(1, 0))   # NA-compatible clique of 4 > v = 2
cat(sprintf("clique larger than v: hers %d constant rows, ours %d\n", nconst(maxconstant(D2)), nconst(maxconstant_c(D2))))
cat(sprintf("\nNA cases %d, NA mismatches %d, cases where her clique exceeds v %d, cases where her rows repeat a constant %d\n", length(nacases), nafails, nabig, nadup))
cat(sprintf("\ncases %d, mismatches %d, total time hers %.1f s, ours %.1f s\n", length(cases) + length(nacases), fails, tot_h, tot_c))
cat("MAXCONSTANT TEST DONE\n")
