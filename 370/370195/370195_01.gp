M=22;

\\ a(n) = Sum_{k=0..floor(n/2)} binomial(2*n,k) * binomial(2*n,n-2*k).
b(n) = sum(k=0, n\2, binomial(2*n, k) * binomial(2*n, n-2*k));
for(n=0, M, print1(b(n), ", "))

a(n, s=2, t=2, u=2) = sum(k=0, n\s, binomial(t*n, k)*binomial(u*n, n-s*k));
for(n=0, M, print1(a(n)-b(n), ", "))

