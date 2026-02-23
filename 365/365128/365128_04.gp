\\ a(n) = (1/(n+1)) * Sum_{k=0..floor(n/3)} (-1)^k * binomial(3*n+3,k) * binomial(4*n-3*k+2,3*n+2).
a(n) = (1/(n+1)) * sum(k=0, n\3, (-1)^k * binomial(3*n+3,k) * binomial(4*n-3*k+2,3*n+2));
for(n=0, 21, print1(a(n), ", "));