M=19;

\\ a(n) = Sum_{k=0..n} binomial(2*n+5*k+1,k) * binomial(2*n+4*k+2,n-k)/(n+2*k+1).
a(n) = sum(k=0, n, binomial(2*n+5*k+1, k) * binomial(2*n+4*k+2, n-k)/(n+2*k+1));
for(n=0, M, print1(a(n),", "))

