M=17;

a(n, r=2, t=4, u=2) = r*sum(k=0, n, binomial(n, k)*binomial(t*n+u*k+r, n)/(t*n+u*k+r));
for(n=0, M, print1(a(n), ", "))

\\ a(n) = Sum_{k=0..n} binomial(n,k) * binomial(4*n+2*k+2,n)/(2*n+k+1).
b(n) = sum(k=0, n, binomial(n, k)*binomial(4*n+2*k+2, n)/(2*n+k+1));
for(n=0, 40, print1(a(n)-b(n),", ")) 


