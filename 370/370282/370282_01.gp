M=19;

\\ a(n) = Sum_{k=0..n} binomial(n+k-1,k) * binomial(4*n+2*k-1,n-k).
a(n) = sum(k=0, n, binomial(n+k-1, k)*binomial(4*n+2*k-1, n-k));
for(n=0, M, print1(a(n), ", "))

