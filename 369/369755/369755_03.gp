bell(n, x) = sum(k=0, n, x^k * stirling(n, k, 2));
\\ a(n) = Sum_{k=0..n} 4^k * Stirling1(n,k) * Bell_k(-1/4)
a(n) = sum(k=0, n, 4^k * stirling(n, k, 1) * bell(k, -1/4));
for(n=0, 20, print1(a(n), ", "));