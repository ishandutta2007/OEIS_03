M=18;

a(n, r=3, t=4, u=1) = r*sum(k=0, n, binomial(n, k)*binomial(t*n+u*k+r, n)/(t*n+u*k+r));
for(n=0, M, print1(a(n),", "))                 

\\ a(n) = 3 * Sum_{k=0..n} binomial(n,k) * binomial(4*n+k+3,n)/(4*n+k+3).
b(n) = 3 * sum(k=0, n, binomial(n,k)*binomial(4*n+k+3,n)/(4*n+k+3));
for(n=0, 100, print1(a(n)-b(n),", "))


