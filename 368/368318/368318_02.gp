M=24;

\\ a(n) = (1/(n+1)) * Sum_{k=0..floor(n/3)} binomial(2*n+2,k) * binomial(2*n+2,n-3*k).
b(n) = (1/(n+1)) * sum(k=0, n\3, binomial(2*n+2,k)*binomial(2*n+2,n-3*k));
for(n=0, M, print1(b(n),", "))  

a(n, s=3, t=2, u=2) = sum(k=0, n\s, binomial(t*(n+1), k)*binomial(u*(n+1), n-s*k))/(n+1);
for(n=0, M, print1(a(n)-b(n),", "))

