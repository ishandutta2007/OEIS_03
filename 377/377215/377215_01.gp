M=28;

\\ a(n) = Sum_{k=0..floor(n/2)} (-4)^k * binomial(-5/2,k) * binomial(n-k-1,n-2*k).
a(n) = sum(k=0, n\2, (-4)^k*binomial(-5/2, k)*binomial(n-k-1, n-2*k));
for(n=0, M, print1(a(n), ", ")) 