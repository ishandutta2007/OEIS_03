M=23;

\\ a(n) = (1/(n+1)) * Sum_{k=0..floor(n/4)} (-1)^k * binomial(n+1,k) * binomial(3*n-3*k+3,n-4*k).
a(n) = sum(k=0, n\4, (-1)^k * binomial(n+1,k) * binomial(3*n-3*k+3,n-4*k)/(n+1));
for(n=0, M, print1(a(n)", "))