M=23;

\\ a(n) = Sum_{k=0..floor(n/2)} binomial(2*n+k-1,k) * binomial(2*n+k,n-2*k).
a(n) = sum(k=0, n\2, binomial(2*n+k-1, k) * binomial(2*n+k, n-2*k));
for(n=0, M, print1(a(n),", "))

\\ a(n) == 0 (mod 2) for n>0.
for(n=1, 100, print1(a(n)%2, ", "))