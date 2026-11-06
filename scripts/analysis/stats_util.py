"""Small statistical helpers shared by the analysis scripts (no SciPy dependency).

hypergeom_upper(k, n, K, N)   P(X >= k) for X hypergeometric: n draws without replacement from a
                              population of N with K successes; computed in log space with lgamma so that
                              populations of tens of thousands are handled without big-integer arithmetic.
bh(pvalues)                   Benjamini-Hochberg adjusted p-values in the input order.
binom_two_sided(k, n, p)      exact two-sided binomial test: sum of the probabilities of all outcomes that are
                              not more likely than the observed one (the definition used by scipy.stats.binomtest).
"""
import math


def _log_comb(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_upper(k, n, K, N):
    if N <= 0 or n <= 0 or K <= 0 or k <= 0:
        return 1.0
    hi = min(n, K)
    if k > hi:
        return 0.0
    denom = _log_comb(N, n)
    logs = [_log_comb(K, x) + _log_comb(N - K, n - x) - denom for x in range(k, hi + 1)]
    m = max(logs)
    return min(1.0, math.exp(m) * sum(math.exp(v - m) for v in logs))


def bh(pvalues):
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    adj = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * n / rank)
        adj[i] = running
    return adj


def binom_two_sided(k, n, p):
    if n <= 0:
        return 1.0
    logpmf = [_log_comb(n, x) + x * math.log(p) + (n - x) * math.log1p(-p) for x in range(n + 1)]
    ref = logpmf[k] + math.log1p(1e-7)          # relative tolerance for ties, as in scipy
    m = max(logpmf)
    total = sum(math.exp(v - m) for v in logpmf if v <= ref)
    return min(1.0, math.exp(m) * total)
