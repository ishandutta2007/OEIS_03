\\ a(n) = Sum_{k=0..floor(n/2)} (-1)^k * binomial(n-k,k)^2.
a(n) = sum(k=0, n\2, (-1)^k * binomial(n-k, k)^2);
for(n=0, 33, print1(a(n), ", "))

