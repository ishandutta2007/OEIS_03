M=21;

\\ a(n) = Sum_{k=0..n} (-4)^k * binomial(-7/2,k) * binomial(n-1,n-k).
a(n) = sum(k=0, n, (-4)^k*binomial(-7/2, k)*binomial(n-1, n-k));
for(n=0, M, print1(a(n), ", ")) 