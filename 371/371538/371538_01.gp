M=20;

\\ a(n) = 3 * Sum_{k=0..n} (-1)^(n-k) * binomial(n-1,n-k) * binomial(6*k+3,k)/(6*k+3).
a(n) = 3 * sum(k=0, n, (-1)^(n-k) * binomial(n-1, n-k) * binomial(6*k+3, k)/(6*k+3));
for(n=0, M, print1(a(n),", "))  

