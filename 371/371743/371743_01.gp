\\ a(n) = Sum_{k=0..floor(n/2)} binomial(4*n-k,n-2*k).
a(n) = sum(k=0, n\2, binomial(4*n-k, n-2*k));
for(n=0, 20, print1(a(n), ", "))

