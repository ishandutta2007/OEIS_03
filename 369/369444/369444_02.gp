M=31;

\\ a(n) = (1/(n+1)) * Sum_{k=0..floor(n/4)} binomial(n+1,k) * binomial(n+1,n-4*k).
b(n) = (1/(n+1)) * sum(k=0, n\4, binomial(n+1,k)*binomial(n+1,n-4*k));
for(n=0, M, print1(b(n),", "))  

a(n, s=4, t=1, u=1) = sum(k=0, n\s, binomial(t*(n+1), k)*binomial(u*(n+1), n-s*k))/(n+1);
for(n=0, M, print1(a(n)-b(n),", "))

