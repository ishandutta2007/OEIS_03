M=18;

\\ a(n) = Sum_{k=0..n} binomial(4*n,k) * binomial(4*n-k-1,n-k).
a(n) = sum(k=0, n, binomial(4*n, k) * binomial(4*n-k-1, n-k));
for(n=0, M, print1(a(n), ", "))

