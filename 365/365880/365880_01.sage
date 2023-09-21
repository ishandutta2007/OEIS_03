def A365880(n):
    h = binomial(6*n + 4, n) * hypergeometric([-n, 4*(n + 1)], [-6 * n - 4], -1) / (n + 1)
    return simplify(h)
print([A365880(n) for n in range(23)])