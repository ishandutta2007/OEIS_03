\\ a(n) = Sum_{k=0..floor(n/3)} (-1)^k * binomial(2*n+2,n-3*k).
a(n) = sum(k=0, n\3, (-1)^k * binomial(2*n+2, n-3*k));
for(n=0, 25, print1(a(n), ", "))

