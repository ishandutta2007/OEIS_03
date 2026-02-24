M=26;

\\ a(n) = (1/(n+1)) * Sum_{k=0..floor(n/2)} (-1)^k * binomial(2*n+k+1,k) * binomial(3*n-2*k+1,n-2*k).
a(n) = (1/(n+1)) * sum(k=0, n\2, (-1)^k * binomial(2*n+k+1,k) * binomial(3*n-2*k+1,n-2*k));
for(n=0, M, print1(a(n),", "));

\\ a(n) = (1/(n+1)) * Sum_{k=0..floor(n/4)} binomial(2*n+k+1,k) * binomial(2*n+2,n-4*k).
b(n) = (1/(n+1)) * sum(k=0, n\4, binomial(2*n+k+1,k) * binomial(2*n+2,n-4*k));
for(n=0, M, print1(a(n)-b(n),", "));



