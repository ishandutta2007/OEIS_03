\\ a(n) = Sum_{k=0..floor(n/2)} binomial(n+k-1,k) * binomial(2*n+k-1,n-2*k).
a(n) = sum(k=0, n\2, binomial(n+k-1,k) * binomial(2*n+k-1,n-2*k));
for(n=0, 23, print1(a(n), ", "))

