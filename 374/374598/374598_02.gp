\\ a(n) = Sum_{k=0..floor(n/2)} 2^k * binomial(n-2*k,k) * binomial(2*(n-2*k),n-2*k).
a(n) = sum(k=0, n\2, 2^k * binomial(n-2*k, k) * binomial(2*(n-2*k), n-2*k));
for(n=0, 50, print1(a(n), ", "))