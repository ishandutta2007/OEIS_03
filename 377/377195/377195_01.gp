M=21;

\\ a(n) = Sum_{k=0..n} (-4)^k * binomial(-5/2,k) * binomial(2*k,n-k).
a(n) = sum(k=0, n, (-4)^k*binomial(-5/2, k)*binomial(2*k, n-k));
for(n=0, M, print1(a(n), ", ")) 