/* pca_core.c - the v rows of a strength 2 array on which the most columns show v distinct symbols
 * (the k1 block of a partitioned covering array, Colbourn et al. 2006). Exact depth-first search over
 * row sets with a column bitset, bounded by the number of columns still alive; stops at a time limit
 * and reports whether it finished. Own code, C. Smolen. */
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

static int N, K, V, W;                 /* rows, columns, symbols, 64-bit words per column set */
static const int *A;                   /* N x K, column-major (R matrix), symbols 0..V-1, NA_INTEGER allowed */
static uint64_t *stackmask;            /* (V+1) x W */
static int *cur, *best, bestk;
static double tlimit; static clock_t t0; static int timedout; static long nodes;

static inline int pop(const uint64_t *m) { int c = 0; for (int w = 0; w < W; w++) c += __builtin_popcountll(m[w]); return c; }

static void dfs(int depth, int from) {
    const uint64_t *m = stackmask + (size_t)depth * W;
    if (depth == V) { int c = pop(m); if (c > bestk) { bestk = c; memcpy(best, cur, V * sizeof(int)); } return; }
    if ((++nodes & 1023) == 0 && (double)(clock() - t0) / CLOCKS_PER_SEC > tlimit) timedout = 1;
    if (timedout) return;
    uint64_t *nm = stackmask + (size_t)(depth + 1) * W;
    for (int r = from; r <= N - (V - depth); r++) {
        memcpy(nm, m, W * sizeof(uint64_t));
        for (int j = 0; j < K; j++) {              /* drop the columns where row r repeats a chosen symbol or is NA */
            if (!(nm[j >> 6] >> (j & 63) & 1)) continue;
            int x = A[(size_t)j * N + r], bad = (x == NA_INTEGER);
            for (int d = 0; d < depth && !bad; d++) if (A[(size_t)j * N + cur[d]] == x) bad = 1;
            if (bad) nm[j >> 6] &= ~((uint64_t)1 << (j & 63));
        }
        if (pop(nm) <= bestk) continue;            /* cannot beat the best row set found */
        cur[depth] = r; dfs(depth + 1, r + 1);
        if (timedout) return;
    }
}

SEXP C_pca_rows(SEXP code, SEXP sv, SEXP sseconds) {
    N = nrows(code); K = ncols(code); V = asInteger(sv); A = INTEGER(code); W = (K + 63) / 64;
    tlimit = asReal(sseconds); t0 = clock(); timedout = 0; nodes = 0; bestk = -1;
    stackmask = (uint64_t *) R_alloc((size_t)(V + 1) * W, sizeof(uint64_t));
    cur = (int *) R_alloc(V, sizeof(int)); best = (int *) R_alloc(V, sizeof(int));
    memset(stackmask, 0, (size_t)W * sizeof(uint64_t));
    for (int j = 0; j < K; j++) stackmask[j >> 6] |= (uint64_t)1 << (j & 63);
    for (int i = 0; i < V; i++) best[i] = i;
    if (N >= V) dfs(0, 0);
    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SEXP rows = PROTECT(allocVector(INTSXP, V));
    for (int i = 0; i < V; i++) INTEGER(rows)[i] = best[i] + 1;
    SET_VECTOR_ELT(out, 0, rows); SET_VECTOR_ELT(out, 1, ScalarInteger(bestk < 0 ? 0 : bestk));
    SET_VECTOR_ELT(out, 2, ScalarLogical(!timedout));
    UNPROTECT(2); return out;
}
