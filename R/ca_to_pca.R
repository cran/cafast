CA_to_PCA_c <- function(D, seconds = 60, ...) {
    ## The partitioned form of a strength 2 array used by CAs::productPCA: v rows on top; in the first
    ## k1 columns they read as the levels in ascending order; in the other k2 columns the level in row i
    ## is not larger than the i-th level. The v rows are chosen to make k1 as large as possible.
    stopifnot(is.matrix(D), is.numeric(D), is.numeric(seconds), length(seconds) == 1, seconds > 0)
    attrs <- attributes(D)
    N <- nrow(D); k <- ncol(D)
    levs <- sort(unique(as.vector(D))); levs <- levs[!is.na(levs)]; v <- length(levs)
    if (N < v) stop("D has fewer rows than levels")
    code <- matrix(match(D, levs) - 1L, N, k); storage.mode(code) <- "integer"
    res <- .Call(C_pca_rows, code, as.integer(v), as.numeric(seconds))
    rows <- res[[1]]; k1 <- res[[2]]; exact <- res[[3]]
    code <- code[c(rows, setdiff(seq_len(N), rows)), , drop = FALSE]
    top <- code[seq_len(v), , drop = FALSE]
    good <- which(apply(top, 2, function(x) !anyNA(x) && length(unique(x)) == v))
    ord <- c(good, setdiff(seq_len(k), good))
    code <- code[, ord, drop = FALSE]
    for (j in seq_len(k)) {
        col <- code[, j]; topj <- col[seq_len(v)]
        if (j <= length(good)) perm <- order(topj) - 1L          # symbol s goes to the row it sits in
        else {                                                 # symbols in order of first appearance: row i gets at most level i
            seen <- unique(topj[!is.na(topj)]); rest <- setdiff(0:(v - 1L), seen)
            perm <- integer(v); perm[c(seen, rest) + 1L] <- 0:(v - 1L)
        }
        code[, j] <- perm[col + 1L]
    }
    out <- matrix(levs[code + 1L], N, k)
    k1 <- length(good); k2 <- k - k1
    type <- if (k1 == 0) NA_character_ else if (k2 == 0 || all(out[seq_len(v), (k1 + 1):k] == levs[1])) "SCA" else "PCA"
    keep <- attrs[setdiff(names(attrs), c("dim", "dimnames", "PCAstatus"))]
    for (a in names(keep)) attr(out, a) <- keep[[a]]
    attr(out, "PCAstatus") <- if (k1 > 0) list(type = type, k1 = k1, k2 = k2) else FALSE
    attr(out, "PCA_rows") <- rows                 # the rows of D that became the top block
    attr(out, "PCA_columns") <- ord               # the columns of D in their new order
    attr(out, "PCA_exact") <- exact               # FALSE when the search stopped at the time limit
    out
}
