M=20;

\\ a(n) = Sum_{k=0..floor(n/3)} binomial(n+k-1,k) * binomial(4*n-1,n-3*k).
a(n) = sum(k=0, n\3, binomial(n+k-1, k)*binomial(4*n-1, n-3*k));
for(n=0, M, print1(a(n), ", "))

