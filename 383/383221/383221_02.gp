M=19;

\\ a(n) = Sum_{k=3..n} 2^(k-3) * 3^(n-k) * binomial(k,3) * |Stirling1(n,k)|.
a(n) = sum(k=3, n, 2^(k-3) * 3^(n-k) * binomial(k, 3) * abs(stirling(n, k, 1)));
for(n=0, M, print1(a(n), ", "));

\\ a(n) = Sum_{k=3..n} (3*n-1)^(k-3) * 3^(n-k) * binomial(k,3) * Stirling1(n,k).
b(n) = sum(k=3, n, (3*n-1)^(k-3) * 3^(n-k) * binomial(k, 3) * stirling(n, k, 1));
for(n=0, 100, print1(a(n)-b(n), ", "));



